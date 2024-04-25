SPHINXOPTS  ?=
SPHINXBUILD ?= sphinx-build
SOURCEDIR   = source
BUILDDIR    = build

.PHONY: Makefile

help:
	@echo "all   : build"
	@echo "clean : clean-up"

all:
	$(SPHINXBUILD) $(SOURCEDIR) $(BUILDDIR)

clean:
	$(RM) -r $(BUILDDIR)
