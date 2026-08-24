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

  - title: Ion Image Visualization
    details: View the spatial intensity distribution of a specified m/z value, and use TIC images and mass spectra to inspect processing results and tissue spatial characteristics.
    link: /ion_image_visualization
    linkText: View Guide


  - title: KMeans Tissue Region Segmentation
    details: Divide pixels into several classes based on low-dimensional features, inspect the resulting spatial regions, and generate removal or retention masks after manual confirmation.
    link: /segmentation
    linkText: View Guide
---
