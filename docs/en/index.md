---
layout: home

hero:
  name: "MassFlow"
  text: "Framework"
  tagline: "A Data Processing Framework for Mass Spectrometry Imaging"
  actions:
    - theme: brand
      text: Getting Started
      link: /getting-started

features:
  - title: Baseline Correction
    details: Remove background signals that change slowly with m/z, making true peaks easier to distinguish from the background and preparing the data for subsequent peak picking and quantitative analysis.
    link: /baseline_correction
    linkText: View Guide

  - title: Noise Reduction
    details: Reduce random fluctuations in mass spectra while preserving true peak shapes as much as possible, improving spectral clarity.
    link: /noise_reduction
    linkText: View Guide

  - title: Intensity Normalization
    details: Reduce differences in overall signal intensity between pixels, making mass-spectrometry signals from different spatial locations easier to compare.
    link: /normalization
    linkText: View Guide

  - title: Peak Picking
    details: Identify major ion peaks from profile spectra and retain only their positions and corresponding intensities, reducing the amount of data for subsequent peak alignment and filtering.
    link: /peak_pick
    linkText: View Guide

  - title: Peak Alignment
    details: Map nearby peaks from different pixels onto a common reference m/z axis so that the same ion can be compared directly across pixels.
    link: /peak_alignment
    linkText: View Guide

  - title: Isotope Peak Filtering
    details: Identify isotope peak groups from the same ion and retain the main peak, reducing repeated contributions from isotope peaks in downstream analysis.
    link: /isotope_filter
    linkText: View Guide

  - title: Adduct Annotation and Filtering
    details: Use mass relationships and spatial distributions to identify different adduct forms that may originate from the same compound, then filter them after manual review.
    link: /adduct_filter
    linkText: View Guide

  - title: Statistical Peak Filtering
    details: Filter low-frequency or low-quality peaks based on their occurrence frequency and signal intensity across pixels, further reducing uninformative features.
    link: /stat_filter
    linkText: View Guide

  - title: Ion Image Visualization
    details: View the spatial intensity distribution of a specified m/z value, and use TIC images and mass spectra to inspect processing results and tissue spatial characteristics.
    link: /ion_image_visualization
    linkText: View Guide

  - title: UMAP Dimensionality Reduction and Visualization
    details: Compress high-dimensional MSI spectral features into a low-dimensional space and use continuous colors to inspect overall spectral differences between pixels.
    link: /umap
    linkText: View Guide

  - title: KMeans Tissue Region Segmentation
    details: Divide pixels into several classes based on low-dimensional features, inspect the resulting spatial regions, and generate removal or retention masks after manual confirmation.
    link: /segmentation
    linkText: View Guide
  
  - title: Best Practices
    details: Introduce the basic Jupyter operations needed before tissue region segmentation, including running cells in order, reusing variables, and rerunning cells after changing parameters.
    link: /jupyter
    linkText: View Guide
  
  - title: R Environment Setup (Optional)
    details: Configure the R environment only if you need to use the R/Cardinal interface. Standard MassFlow preprocessing does not depend on R, so most users can skip this page.
    link: /r_environment_setup
    linkText: View Guide
---
