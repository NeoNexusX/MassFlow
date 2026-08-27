import { defineConfig, type HeadConfig } from 'vitepress'
import MarkdownItMathJaX3PRO from 'markdown-it-mathjax3-pro'
import { withMermaid } from 'vitepress-plugin-mermaid'

/**
 * Copyright year range shown in the footer.
 *
 * Resolved from the system clock when the site is built, so the footer no
 * longer needs to be edited by hand each year. Collapses to a single year
 * while the current year is still the starting year, avoiding '2025-2025'.
 *
 * Parameters: None.
 * Returns: Either '<start>' or '<start>-<current>'.
 * Exceptions: None.
 */
const COPYRIGHT_START_YEAR = 2025
const currentYear = new Date().getFullYear()
const copyrightYears =
  currentYear > COPYRIGHT_START_YEAR
    ? `${COPYRIGHT_START_YEAR}-${currentYear}`
    : `${COPYRIGHT_START_YEAR}`

/**
 * Site Config with i18n locales
 *
 * This configuration enables bilingual documentation with English as root and Chinese under /zh.
 * To satisfy the current VitePress version typings, per-locale theme configuration
 * is not used; instead, a single themeConfig is defined.
 *
 * Exceptions: None.
 */
export default withMermaid(
  defineConfig({
  /**
   * Base path for GitHub Pages deployment
   *
   * Parameters: None.
   * Returns: Static configuration object field.
   * Exceptions: None.
   */
  base: '/MassFlow/',
  title: 'MassFlow',
  description: 'MassFlow',
  cleanUrls: true,
  rewrites: {
    /** Map English content folder to root paths */
    'en/:path': ':path'
  },
  locales: {
    /** English as the root locale */
    root: {
      label: 'English',
      lang: 'en-US',

      themeConfig: {
          nav: [
            { text: 'Home', link: '/' },
            { text: 'Getting Started', link: '/getting-started' }
          ]
        }
    },
    /** Simplified Chinese locale under /zh/ */
    zh: {
      label: '简体中文',
      lang: 'zh-CN',

      themeConfig: {
          nav: [
            { text: '首页', link: '/zh/' },
            { text: '快速开始', link: '/zh/getting-started' }
          ]
        }
    }
  },
  mermaid: {
      // refer https://mermaid.js.org/config/setup/modules/mermaidAPI.html#mermaidapi-configuration-defaults for options
  },
    // optionally set additional config for plugin itself with MermaidPluginConfig
  mermaidPlugin: {
      class: "mermaid", // set additional css classes for parent container
  },
  markdown: {
    /**
     * Enable built-in Mermaid support for fenced code blocks.
     * 
     * Parameters: None.
     * Returns: Static configuration property enabling Mermaid rendering.
     * Exceptions: None. If Mermaid library is missing, VitePress will log an error during dev/build.
     */
    config: (md) => {
      md.use(MarkdownItMathJaX3PRO, {
        // add new inlineMathSeparator && displayMathSeparator
        tex: {
          inlineMath: [['$', '$'], ['§', '§']],
          displayMath: [['$$', '$$'], ['§§', '§§']],
        },
        //enable chtml mode 
        chtml: {
          fontURL: 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/output/chtml/fonts/woff-v2'
        }
      })

    },
    theme: {
      light: 'catppuccin-latte',
      dark: 'catppuccin-macchiato'
    }
  },

  transformPageData(pageData) {
    const head = (pageData.frontmatter.head ??= []);
    const inject_content = pageData.frontmatter.inject_content;
    if (inject_content && Array.isArray(inject_content)) {
      inject_content.forEach(item => {
        const { type, contribution, content } = item;
        const headEntry = [type, contribution || {}, content || ''].filter(Boolean);
        head.push(headEntry as HeadConfig);
      });
      delete pageData.frontmatter.inject_content;
    }
  },

  themeConfig: {
    // Single theme config (compatible with current typings)

    /**
     * Sidebar configuration
     *
     * This lists all currently available documentation pages in both English and Chinese locales.
     *
     * Parameters: None.
     * Returns: Not applicable (static configuration object).
     * Exceptions: None.
     */
    sidebar: {
  '/zh/': [
    {
      text: '简体中文',
      items: [
        { text: '首页', link: '/zh/' },
        { text: '快速开始', link: '/zh/getting-started' },
        { text: '最佳实践', link: '/zh/jupyter'},
        { text: '预处理方法',
          collapsed: false,
          items:[
        { text: '基线校正', link: '/zh/baseline_correction' },
        { text: '噪声降低', link: '/zh/noise_reduction' },
        { text: '强度归一化', link: '/zh/normalization' },
        { text: '峰提取', link: '/zh/peak_pick' },
        { text: '峰对齐', link: '/zh/peak_alignment' },
        { text: '同位素峰过滤', link: '/zh/isotope_filter' },
        { text: '加合物注释与过滤', link: '/zh/adduct_filter' },
        { text: '统计峰过滤', link: '/zh/stat_filter' }
          ]
        },
        {
          text: '空间分析',
          collapsed: false,
          items: [
            { text: '离子图可视化', link: '/zh/ion_image_visualization'},
            { text: 'UMAP', link: '/zh/umap'},
            { text: 'KMeans 聚类与区域选择', link: '/zh/segmentation'},
          ]
        },
        {
          text: '环境配置',
          collapsed: false,
          items: [
            {text: 'R 环境配置', link: '/zh/r_environment_setup'}
          ]
        }
      ]
    }
  ],

  '/': [
    {
      text: 'English',
      items: [
        { text: 'Home', link: '/' },
        { text: 'Getting Started', link: '/getting-started' },
        { text: 'Best Practices', link: '/jupyter' },
        {
          text: 'Preprocessing Methods',
          collapsed: false,
          items: [
            {text: 'Baseline Correction', link: '/baseline_correction'},
            {text: 'Noise Reduction', link: '/noise_reduction'},
            {text: 'Intensity Normalization', link: '/normalization'},
            {text: 'Peak Picking', link: '/peak_pick'},
            {text: 'Peak Alignment', link: '/peak_alignment'},
            {text: 'Isotope Peak Filtering', link: '/isotope_filter'},
            {text: 'Adduct Annotation and Filtering', link: '/adduct_filter'},
            {text: 'Statistical Peak Filtering', link: '/stat_filter'},
          ]
        },
          {
          text: 'Spatial Analysis',
            collapsed: false,
          items: [
        { text: 'Ion Image Visualization', link: '/ion_image_visualization' },
        { text: 'UMAP', link: '/umap' },
        { text: 'KMeans Clustering and Region Selection', link: '/segmentation' },

              ]
          },
        {
          text: 'Environment Setup',
          collapsed: false,
          items: [
            {text: 'R Environment Setup', link: '/r_environment_setup'}
          ]
        }
      ]
    }
  ]
},
    socialLinks: [
      { icon: 'github', link: 'https://github.com/NeoNexusX/MassFlow' }
    ],
    footer: {
      message: 'Released under the GNU License.',
      copyright: `Copyright © ${copyrightYears} <a href="https://bionet.xmu.edu.cn/">Bionet Team</a>`
    },
    search: {
      provider: 'local'
    }
  }
})
)

