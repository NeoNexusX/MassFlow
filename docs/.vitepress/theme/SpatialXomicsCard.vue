<script setup lang="ts">
/**
 * Aside promo card linking the docs to SpatialXomics, the open spatial-omics
 * data platform maintained by the same team.
 *
 * The card ships folded away into its own viridis strip: all that shows in
 * the aside is a slowly drifting gradient slit, which doubles as the control
 * that opens the card. The choice is remembered in localStorage, so a reader
 * who opened or dismissed it keeps that state on the next page.
 *
 * Every colour is a VitePress theme variable so the card follows the
 * light/dark palette automatically -- the one exception is the viridis strip,
 * which is a scientific colormap rather than a theme colour and therefore
 * stays fixed in both themes.
 */
import { onMounted, ref } from 'vue'

const HREF = 'https://sxo.bionet.pro/'
const STORAGE_KEY = 'massflow-sxo-card-collapsed'

/**
 * Starts folded: the card is an aside promo, so it opens only for a reader
 * who asked for it, and the choice then sticks.
 */
const collapsed = ref(true)

/**
 * Gates the open/close transition until after the stored state has been
 * applied, so a card left collapsed appears already folded instead of
 * visibly folding itself away on every page load.
 */
const animated = ref(false)

function toggle() {
  collapsed.value = !collapsed.value
  localStorage.setItem(STORAGE_KEY, String(collapsed.value))
}

onMounted(() => {
  collapsed.value = localStorage.getItem(STORAGE_KEY) !== 'false'
  requestAnimationFrame(() => {
    animated.value = true
  })
})
</script>

<template>
  <section
    class="sxo"
    :class="{ 'is-collapsed': collapsed, 'is-animated': animated }"
  >
    <button
      class="sxo-spectrum"
      type="button"
      :aria-expanded="!collapsed"
      aria-controls="sxo-panel"
      :aria-label="
        collapsed
          ? 'Expand the SpatialXomics card'
          : 'Collapse the SpatialXomics card'
      "
      :title="collapsed ? 'Expand' : 'Collapse'"
      @click="toggle"
    >
      <svg
        class="sxo-chevron"
        viewBox="0 0 12 7"
        aria-hidden="true"
        focusable="false"
      >
        <path
          d="M1.5 5 6 1.6 10.5 5"
          fill="none"
          stroke="currentColor"
          stroke-width="2"
          stroke-linecap="round"
          stroke-linejoin="round"
        />
      </svg>
    </button>

    <!-- inert keeps the clipped link out of the tab order while folded. -->
    <div id="sxo-panel" class="sxo-panel" :inert="collapsed">
      <div class="sxo-panel-inner">
        <a
          class="sxo-link"
          :href="HREF"
          target="_blank"
          rel="noreferrer"
          aria-label="SpatialXomics, an open platform for spatial omics data (opens in a new tab)"
        >
          <span class="sxo-eyebrow">Related project</span>

          <span class="sxo-wordmark">Spatial<span class="sxo-x">X</span>omics</span>

          <span class="sxo-desc">
            An open platform for spatial omics data.
          </span>
          <span class="sxo-desc">
            Browse, download, and explore
            datasets. Free for everyone.
          </span>

          <span class="sxo-cta">
            View datasets<span class="sxo-arrow" aria-hidden="true">&rarr;</span>
          </span>
        </a>
      </div>
    </div>
  </section>
</template>

<style scoped>
.sxo {
  display: block;
  overflow: hidden;
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  background-color: var(--vp-c-bg-elv);
  color: inherit;
  transition:
    border-color 0.25s,
    border-radius 0.25s;
}

.sxo:hover {
  border-color: var(--vp-c-brand-1);
}

/* A ten-pixel slit reads better with a tighter corner than the full card. */
.sxo.is-collapsed {
  border-radius: 6px;
}

/*
 * Signature element, and the fold control: the viridis colormap, shared by the
 * ion images in these docs and by the SpatialXomics landing page. Fixed stops
 * on purpose -- a colormap that shifted with the theme would no longer be
 * viridis.
 */
.sxo-spectrum {
  position: relative;
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 0 12px;
  width: 100%;
  height: 15px;
  border: 0;
  background: none;
  overflow: hidden;
  cursor: pointer;
}

/*
 * The drifting colormap lives on a track twice the strip's width, carrying two
 * identical tiles. Each tile runs viridis out to yellow and back again, so its
 * two edges are the same purple and the seam between tiles is invisible;
 * sliding the track by exactly one tile therefore loops without a jump.
 */
.sxo-spectrum::before {
  content: '';
  position: absolute;
  inset: 0;
  width: 200%;
  background-image: linear-gradient(
    90deg,
    #440154,
    #3b528b,
    #21918c,
    #5ec962,
    #fde725,
    #5ec962,
    #21918c,
    #3b528b,
    #440154
  );
  background-size: 50% 100%;
  background-repeat: repeat-x;
  animation: sxo-drift 12s linear infinite;
}

@keyframes sxo-drift {
  from {
    transform: translateX(0);
  }

  to {
    transform: translateX(-50%);
  }
}

.sxo-spectrum:focus-visible {
  outline: 2px solid var(--vp-c-brand-1);
  outline-offset: 2px;
}

/*
 * White, because viridis spends most of its range dark: purple through teal
 * carry the glyph easily, and only the short yellow-green stretch weakens it.
 * No plate behind it -- the drifting colormap stays the whole of the strip.
 */
.sxo-chevron {
  position: relative;
  width: 14px;
  height: 8px;
  color: rgba(255, 255, 255, 0.92);
  transition:
    transform 0.25s,
    color 0.25s;
}

.sxo-spectrum:hover .sxo-chevron {
  color: #fff;
}

.sxo.is-collapsed .sxo-chevron {
  transform: rotate(180deg);
}

/*
 * 1fr -> 0fr collapses the row to nothing without hard-coding the card's
 * height, which depends on how the copy wraps in the aside.
 */
.sxo-panel {
  display: grid;
  grid-template-rows: 1fr;
}

.sxo.is-animated .sxo-panel {
  transition: grid-template-rows 0.3s ease;
}

.sxo.is-collapsed .sxo-panel {
  grid-template-rows: 0fr;
}

.sxo-panel-inner {
  overflow: hidden;
  min-height: 0;
}

.sxo-link {
  display: block;
  padding: 24px;
  color: inherit;
  text-decoration: none;
}

.sxo-eyebrow {
  display: block;
  font-family: var(--vp-font-family-mono);
  font-size: 15px;
  font-weight: 500;
  letter-spacing: 0.12em;
  line-height: 1;
  text-transform: uppercase;
  color: var(--vp-c-text-3);
}

.sxo-wordmark {
  display: block;
  margin-top: 12px;
  font-size: 24px;
  font-weight: 600;
  letter-spacing: -0.01em;
  line-height: 1.2;
  color: var(--vp-c-text-1);
}

/* Mirrors the tinted X of the SpatialXomics wordmark. */
.sxo-x {
  font-weight: 700;
  color: var(--vp-c-brand-1);
}

.sxo-desc {
  display: block;
  margin-top: 10px;
  font-size: 16px;
  line-height: 1.55;
  color: var(--vp-c-text-2);
}

.sxo-cta {
  display: block;
  margin-top: 16px;
  font-size: 16px;
  font-weight: 600;
  line-height: 1;
  color: var(--vp-c-brand-1);
}

.sxo-arrow {
  display: inline-block;
  margin-left: 6px;
  transition: transform 0.25s;
}

.sxo-link:hover .sxo-arrow {
  transform: translateX(3px);
}

@media (prefers-reduced-motion: reduce) {
  .sxo,
  .sxo-chevron,
  .sxo-arrow,
  .sxo.is-animated .sxo-panel {
    transition: none;
  }

  .sxo-spectrum::before {
    animation: none;
  }

  .sxo-link:hover .sxo-arrow {
    transform: none;
  }
}
</style>
