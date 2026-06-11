from __future__ import annotations

import math
import os
from concurrent.futures import ThreadPoolExecutor
from typing import Literal, Optional, Sequence

import numpy as np
from scipy.special import expit
from scipy.sparse import bmat, coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

EPS = 1e-12
WeakComponentFallback = Literal["solve", "keep_zero"]


class RLNormError(ValueError):
    """Raised when input data or parameters are invalid for RL-Norm."""


def _validate_input(mat_data: np.ndarray) -> np.ndarray:
    arr = np.asarray(mat_data, dtype=float)
    if arr.ndim != 3:
        raise RLNormError("mat_data must be a 3D array with shape (n_row, n_col, n_ion)")
    n_row, n_col, n_ion = arr.shape
    if n_row <= 0 or n_col <= 0 or n_ion <= 0:
        raise RLNormError("mat_data dimensions must all be > 0")
    if not np.all(np.isfinite(arr)):
        raise RLNormError("mat_data contains NaN or Inf")
    return np.maximum(arr, 0.0)


def _infer_mask(mat_data: np.ndarray) -> np.ndarray:
    return np.any(mat_data > 0, axis=2)


def _flatten_valid_pixels(mask: np.ndarray) -> tuple[np.ndarray, dict[tuple[int, int], int]]:
    coords = np.argwhere(mask)
    if coords.size == 0:
        raise RLNormError("No valid pixels in mask")
    index_map = {tuple(coord): idx for idx, coord in enumerate(coords)}
    return coords, index_map


def _neighbor_offsets(window_size: int) -> list[tuple[int, int]]:
    radius = window_size // 2
    offsets: list[tuple[int, int]] = []
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            if dr == 0 and dc == 0:
                continue
            offsets.append((dr, dc))
    return offsets


def _generate_candidate_edges(
    coords: np.ndarray,
    index_map: dict[tuple[int, int], int],
    n_row: int,
    n_col: int,
    window_size: int,
) -> list[tuple[int, int, int, int]]:
    offsets = _neighbor_offsets(window_size)
    edges: list[tuple[int, int, int, int]] = []
    for i, (r, c) in enumerate(coords):
        r_int = int(r)
        c_int = int(c)
        for dr, dc in offsets:
            rr = r_int + dr
            cc = c_int + dc
            if rr < 0 or rr >= n_row or cc < 0 or cc >= n_col:
                continue
            j = index_map.get((rr, cc))
            if j is None or j <= i:
                continue
            edges.append((i, j, dr, dc))
    return edges


def _edge_chunk_worker(
    edges_chunk: Sequence[tuple[int, int, int, int]],
    spectra: np.ndarray,
    spectra_nz: np.ndarray,
    spectra_norm: np.ndarray,
    sigma_space: float,
    tau: float,
    gamma: float,
    min_common: int,
    min_weight: float,
    trim_fraction: float,
) -> list[tuple[int, int, float, float]]:
    out: list[tuple[int, int, float, float]] = []
    sigma_denom = 2.0 * sigma_space * sigma_space + EPS

    for i, j, dr, dc in edges_chunk:
        xi = spectra[i]
        xj = spectra[j]

        norm_prod = spectra_norm[i] * spectra_norm[j]
        if norm_prod <= EPS:
            continue

        dist2 = float(dr * dr + dc * dc)
        w_space = math.exp(-dist2 / sigma_denom)
        dot_ij = float(np.dot(xi, xj))
        w_spec = max(float(dot_ij / (norm_prod + EPS)), 0.0)

        nz_i = spectra_nz[i]
        nz_j = spectra_nz[j]
        common = nz_i & nz_j
        union = nz_i | nz_j

        union_count = int(np.sum(union))
        if union_count == 0:
            continue

        common_count = int(np.sum(common))
        if common_count < min_common:
            continue

        jac = common_count / (union_count + EPS)
        penalty = float(expit(gamma * (jac - tau)))
        weight = w_space * w_spec * penalty
        if weight < min_weight:
            continue

        common_xi = xi[common]
        common_xj = xj[common]
        intensity = 0.5 * (common_xi + common_xj)
        order = np.argsort(intensity)

        trim_k = 0
        if len(order) >= 10:
            trim_k = max(1, int(len(order) * trim_fraction))
        if trim_k > 0 and len(order) > 2 * trim_k:
            keep = order[trim_k:-trim_k]
        else:
            keep = order
        if len(keep) == 0:
            continue

        log_diff = np.log(np.maximum(common_xi[keep], EPS)) - np.log(np.maximum(common_xj[keep], EPS))
        rij = float(np.median(log_diff))
        out.append((i, j, weight, rij))

    return out


def _build_rl_graph(
    spectra: np.ndarray,
    edges: Sequence[tuple[int, int, int, int]],
    sigma_space: float,
    tau: float,
    gamma: float,
    min_common: int,
    min_weight: float,
    trim_fraction: float,
    n_jobs: int,
    chunk_size: int,
) -> list[tuple[int, int, float, float]]:
    if len(edges) == 0:
        return []

    spectra_nz = spectra > 0
    spectra_norm = np.linalg.norm(spectra, axis=1)

    if n_jobs <= 1 or len(edges) <= chunk_size:
        return _edge_chunk_worker(
            edges,
            spectra,
            spectra_nz,
            spectra_norm,
            sigma_space,
            tau,
            gamma,
            min_common,
            min_weight,
            trim_fraction,
        )

    chunks: list[Sequence[tuple[int, int, int, int]]] = []
    for start in range(0, len(edges), chunk_size):
        chunks.append(edges[start:start + chunk_size])

    edge_list: list[tuple[int, int, float, float]] = []
    max_workers = min(n_jobs, len(chunks))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _edge_chunk_worker,
                chunk,
                spectra,
                spectra_nz,
                spectra_norm,
                sigma_space,
                tau,
                gamma,
                min_common,
                min_weight,
                trim_fraction,
            )
            for chunk in chunks
        ]
        for fut in futures:
            edge_list.extend(fut.result())

    return edge_list


def _edge_list_to_adjacency(n_pixels: int, edge_list: Sequence[tuple[int, int, float, float]]) -> csr_matrix:
    if len(edge_list) == 0:
        return csr_matrix((n_pixels, n_pixels), dtype=float)
    rows = np.array([e[0] for e in edge_list], dtype=int)
    cols = np.array([e[1] for e in edge_list], dtype=int)
    weights = np.array([e[2] for e in edge_list], dtype=float)
    return csr_matrix(coo_matrix(
        (np.r_[weights, weights], (np.r_[rows, cols], np.r_[cols, rows])),
        shape=(n_pixels, n_pixels),
    ))


def _solve_componentwise_zero_mean(
    n_pixels: int,
    edge_list: Sequence[tuple[int, int, float, float]],
    *,
    weak_component_fallback: WeakComponentFallback,
) -> np.ndarray:
    if len(edge_list) == 0:
        return np.zeros(n_pixels, dtype=float)

    rows = np.array([e[0] for e in edge_list], dtype=int)
    cols = np.array([e[1] for e in edge_list], dtype=int)
    w = np.array([e[2] for e in edge_list], dtype=float)
    r = np.array([e[3] for e in edge_list], dtype=float)

    adj = _edge_list_to_adjacency(n_pixels, edge_list)
    n_comp, labels = connected_components(adj, directed=False)
    a = np.zeros(n_pixels, dtype=float)

    for comp_id in range(n_comp):
        idx = np.where(labels == comp_id)[0]
        m = len(idx)
        if m <= 1:
            continue

        pos = {old: new for new, old in enumerate(idx)}
        sel = np.isin(rows, idx) & np.isin(cols, idx)
        rr = rows[sel]
        cc = cols[sel]
        ww = w[sel]
        rv = r[sel]
        n_edges_comp = len(rr)

        if n_edges_comp < m - 1 and weak_component_fallback == "keep_zero":
            continue

        lr = np.array([pos[int(v)] for v in rr], dtype=int)
        lc = np.array([pos[int(v)] for v in cc], dtype=int)

        data = np.r_[ww, ww, -ww, -ww]
        row_idx = np.r_[lr, lc, lr, lc]
        col_idx = np.r_[lr, lc, lc, lr]
        l_matrix = csr_matrix(coo_matrix((data, (row_idx, col_idx)), shape=(m, m)))

        b = np.zeros(m, dtype=float)
        np.add.at(b, lr, ww * rv)
        np.add.at(b, lc, -ww * rv)

        ones = np.ones((m, 1), dtype=float)
        kkt = bmat([[l_matrix, csr_matrix(ones)], [csr_matrix(ones.T), None]], format="csr")
        rhs = np.r_[-b, 0.0]
        sol = spsolve(kkt, rhs)
        a[idx] = np.asarray(sol[:m]).reshape(-1)

    return a


def rl_norm_cube(
    mat_data: np.ndarray,
    *,
    mask: Optional[np.ndarray] = None,
    window_size: int = 5,
    sigma_space: float = 1.6,
    tau: float = 0.30,
    gamma: float = 8.0,
    trim_fraction: float = 0.10,
    min_common: int = 10,
    min_weight: float = 0.06,
    max_iter: int = 1,
    tol: float = 1e-4,
    n_jobs: Optional[int] = None,
    chunk_size: int = 4000,
    weak_component_fallback: WeakComponentFallback = "solve",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run RL-Norm on one MSI cube and return (normalized_cube, Aij, logA)."""
    if window_size < 3 or window_size % 2 == 0:
        raise RLNormError("window_size must be an odd integer >= 3")
    if sigma_space <= 0:
        raise RLNormError("sigma_space must be > 0")
    if not (0 <= trim_fraction < 0.5):
        raise RLNormError("trim_fraction must satisfy 0 <= trim_fraction < 0.5")
    if min_common < 1:
        raise RLNormError("min_common must be >= 1")
    if min_weight < 0:
        raise RLNormError("min_weight must be >= 0")
    if max_iter < 1:
        raise RLNormError("max_iter must be >= 1")
    if tol < 0:
        raise RLNormError("tol must be >= 0")
    if chunk_size < 1:
        raise RLNormError("chunk_size must be >= 1")
    if weak_component_fallback not in {"solve", "keep_zero"}:
        raise RLNormError("weak_component_fallback must be 'solve' or 'keep_zero'")

    arr = _validate_input(mat_data)
    n_row, n_col, n_ion = arr.shape

    if mask is None:
        mask_arr = _infer_mask(arr)
    else:
        mask_arr = np.asarray(mask, dtype=bool)
        if mask_arr.shape != (n_row, n_col):
            raise RLNormError("mask shape must be (n_row, n_col)")

    n_valid = int(np.sum(mask_arr))
    if n_valid == 0:
        raise RLNormError("mask contains no valid pixel")

    min_common = min(int(min_common), int(n_ion))

    if n_jobs is None:
        n_jobs = max(1, os.cpu_count() or 1)
    else:
        n_jobs = max(1, int(n_jobs))

    coords, index_map = _flatten_valid_pixels(mask_arr)
    candidate_edges = _generate_candidate_edges(coords, index_map, n_row, n_col, window_size)

    loga_valid_total = np.zeros(n_valid, dtype=float)
    mat_iter = arr.copy()

    for _ in range(max_iter):
        spectra = mat_iter[mask_arr]

        edge_list = _build_rl_graph(
            spectra=spectra,
            edges=candidate_edges,
            sigma_space=sigma_space,
            tau=tau,
            gamma=gamma,
            min_common=min_common,
            min_weight=min_weight,
            trim_fraction=trim_fraction,
            n_jobs=n_jobs,
            chunk_size=chunk_size,
        )

        delta_loga_valid = _solve_componentwise_zero_mean(
            n_valid,
            edge_list,
            weak_component_fallback=weak_component_fallback,
        )
        loga_valid_total += delta_loga_valid

        loga = np.zeros((n_row, n_col), dtype=float)
        loga[mask_arr] = loga_valid_total
        aij = np.ones((n_row, n_col), dtype=float)
        aij[mask_arr] = np.exp(loga[mask_arr])
        mat_iter = arr * aij[:, :, None]

        delta_max = float(np.max(np.abs(delta_loga_valid))) if len(delta_loga_valid) > 0 else 0.0
        if delta_max <= tol:
            break

    loga = np.zeros((n_row, n_col), dtype=float)
    loga[mask_arr] = loga_valid_total
    aij = np.ones((n_row, n_col), dtype=float)
    aij[mask_arr] = np.exp(loga[mask_arr])
    mat_y = arr * aij[:, :, None]
    return mat_y, aij, loga
