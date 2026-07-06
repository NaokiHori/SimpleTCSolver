# Documentation

Execute the following command to build a set of HTML files:

```bash
docker run --rm -e PYTHONDONTWRITEBYTECODE=1 -v ${PWD}:/docs sphinxdoc/sphinx sphinx-build source html
```

Open `html/index.html` using your favorite browser.
