# LocusSnap website

This directory contains the static, multi-page project website for GitHub
Pages. The information architecture separates onboarding, task guides,
copyable recipes, configuration, technical reference, the curated gallery,
and troubleshooting:

```text
index.html              project landing page
getting-started.html    installation and first snapshot
workflows.html          task-oriented analysis guides
recipes.html            copy-paste command patterns
configuration.html      YAML, styling, tracks, and plugins
reference.html          CLI, formats, output, and Python API
gallery.html            maintained representative figures
faq.html                common questions and troubleshooting
assets/                 shared CSS, JavaScript, icons, and figures
```

## Preview locally

From the repository root:

```bash
python3 -m http.server 8000 --directory docs
```

Then open <http://localhost:8000/>.

## Publish

In the GitHub repository, open **Settings → Pages**, select **Deploy from a
branch**, and choose the `main` branch with the `/docs` folder. Future pushes to
`main` will update the site automatically.

The pages deliberately share `assets/site.css` and `assets/site.js` so the
landing page and reference guides behave consistently. Gallery images are
copied into `docs/assets/images` because GitHub Pages only publishes files
inside the selected `/docs` source folder.
