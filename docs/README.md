# LocusSnap website

This directory contains the static project website for GitHub Pages.

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

The gallery images are intentionally copied into `docs/assets/images` because
GitHub Pages only publishes files inside the selected `/docs` source folder.
