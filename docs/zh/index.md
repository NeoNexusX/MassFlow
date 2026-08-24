---
layout: home

hero:
  name: "MassFlow"
  text: "框架"
  tagline: "面向质谱成像的数据处理框架"
  actions:
    - theme: brand
      text: 快速开始
      link: /zh/getting-started

features:
  - title: 基线校正
    details: 去除随 m/z 缓慢变化的背景信号，使真实峰更加容易与背景区分，并为后续峰提取和定量分析做好准备。
    link: /zh/baseline_correction
    linkText: 查看指南

  - title: 噪声抑制
    details: 减少质谱中的随机波动，在尽量保留真实峰形的同时提高谱图的清晰度。
    link: /zh/noise_reduction
    linkText: 查看指南

  - title: 强度归一化
    details: 减小不同像素之间整体信号强弱造成的差异，便于比较不同空间位置的质谱信号。
    link: /zh/normalization
    linkText: 查看指南

  - title: 峰提取
    details: 从轮廓谱中识别主要离子峰，只保留峰位置及其对应强度，为后续峰对齐和峰过滤减少数据量。
    link: /zh/peak_pick
    linkText: 查看指南

  - title: 峰对齐
    details: 将不同像素中位置接近的峰对应到统一的参考 m/z 轴，使同一种离子能够在不同像素之间直接比较。
    link: /zh/peak_alignment
    linkText: 查看指南

  - title: 同位素峰过滤
    details: 识别同一离子的同位素峰组，并保留主要峰，减少后续同位素峰对分析结果的重复影响。
    link: /zh/isotope_filter
    linkText: 查看指南

  - title: 加合物注释与过滤
    details: 根据质量关系和空间分布识别可能来自同一化合物的不同加合物形式，并在人工检查后进行过滤。
    link: /zh/adduct_filter
    linkText: 查看指南

  - title: 统计峰过滤
    details: 根据峰在不同像素中的出现频率和信号强度筛选低频或低质量峰，进一步减少无效特征。
    link: /zh/stat_filter
    linkText: 查看指南

  - title: 离子图可视化
    details: 查看指定 m/z 对应的空间强度分布，并结合 TIC 和质谱图检查处理结果及组织空间特征。
    link: /zh/ion_image_visualization
    linkText: 查看指南

  - title: UMAP 降维与可视化
    details: 将高维 MSI 质谱特征压缩到低维空间，用连续颜色观察像素之间整体质谱差异。
    link: /zh/umap
    linkText: 查看指南

  - title: KMeans 组织区域划分
    details: 基于低维特征将像素划分为若干类别，查看空间区域效果，并在人工确认后生成删除或保留掩码。
    link: /zh/segmentation
    linkText: 查看指南
  
  - title: 最佳实践
    details: 介绍组织区域划分前需要掌握的 Jupyter 基本操作，包括按顺序运行代码单元格、复用变量和修改参数后重新运行。
    link: /zh/jupyter
    linkText: 查看指南
  
  - title:  R 环境设置（可选）
    details: 仅在需要使用 R/Cardinal 接口时配置 R 环境。普通 MassFlow 预处理不依赖 R，可直接跳过本页。
    link: /zh/r_environment_setup
    linkText: 查看指南
---