import { defineConfig, type HeadConfig } from 'vitepress'
import MarkdownItMathJaX3PRO from 'markdown-it-mathjax3-pro'
import { withMermaid } from 'vitepress-plugin-mermaid'
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
      lang: 'en-US'
    },
    /** Simplified Chinese locale under /zh/ */
    zh: {
      label: '简体中文',
      lang: 'zh-CN'
    }
  },
  mermaid: {
      // refer https://mermaid.js.org/config/setup/modules/mermaidAPI.html#mermaidapi-configuration-defaults for options
  },
    // optionally set additional config for plugin itself with MermaidPluginConfig
  mermaidPlugin: {
      class: "mermaid my-class", // set additional css classes for parent container 
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
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/getting-started' },
      { text: 'Contribution', link: '/contribution' }
    ],
    /**
     * Sidebar configuration
     *
     * This lists all currently available documentation pages in both English and Chinese locales.
     *
     * Parameters: None.
     * Returns: Not applicable (static configuration object).
     * Exceptions: None.
     */
    sidebar: [
      {
        text: 'English',
        items: [
          { text: 'Home', link: '/' },
          { text: 'Getting Started', link: '/getting-started' },
          { text: 'Contribution', link: '/contribution' },
          { text: 'Collaboration Guide', link: '/collaboration_guide' },
          { text: 'Naming Conventions', link: '/naming-conventions' },
          { text: 'Data Structures', link: '/ms-data-structures' }
        ]
      },
    ],
    socialLinks: [
      { icon: 'github', link: 'https://github.com/NeoNexusX/MassFlow' }
    ],
    footer: {
      message: 'Released under the GNU License.',
      copyright: 'Copyright © 2025-present <a href="https://bionet.xmu.edu.cn/">Bionet Team</a>'
    },
    search: {
      provider: 'local'
    }
  }
})
)

