---
layout: home

hero:
  name: "MassFlow"
  text: "计算引擎"
  tagline: "面向质谱成像的数据处理引擎"
  actions:
    - theme: brand
      text: 快速开始
      link: /zh/getting-started
  image: https://official-oss.oss-cn-hongkong.aliyuncs.com/docs/dataset_2_umap%20(1).png

features:

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

  - title: 离子图可视化
    details: 查看指定 m/z 对应的空间强度分布，并结合 TIC 和质谱图检查处理结果及组织空间特征。
    link: /zh/ion_image_visualization
    linkText: 查看指南


  - title: KMeans 组织区域划分
    details: 基于低维特征将像素划分为若干类别，查看空间区域效果，并在人工确认后生成删除或保留掩码。
    link: /zh/segmentation
    linkText: 查看指南


---