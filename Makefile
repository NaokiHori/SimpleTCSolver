CC        := mpicc
CFLAGS    := -std=c99 -Wall -Wextra -O3
DEPEND    := -MMD
LIBS      := -lfftw3 -lm
INCLUDES  := -Iinclude
SRCSDIR   := src
OBJSDIR   := obj
SRCS      := $(foreach dir, $(shell find $(SRCSDIR) -type d), $(wildcard $(dir)/*.c))
OBJS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.o)))
DEPS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.d)))
OUTPUTDIR := output
TARGET    := a.out

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files \"$(OBJSDIR)/*.o\""
	@echo "output  : create directories to store artifacts"
	@echo "datadel : remove contents in sub-directories of \"$(OUTPUTDIR)\""
	@echo "help    : show this message"

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(DEPEND) -o $@ $^ $(LIBS)

$(OBJSDIR)/%.o: $(SRCSDIR)/%.c
	@if [ ! -e `dirname $@` ]; then \
		mkdir -p `dirname $@`; \
	fi
	$(CC) $(CFLAGS) $(DEPEND) $(INCLUDES) -c $< -o $@

clean:
	$(RM) -r $(OBJSDIR) $(TARGET)

output:
	@if [ ! -e $(OUTPUTDIR)/log ]; then \
		mkdir -p $(OUTPUTDIR)/log; \
	fi
	@if [ ! -e $(OUTPUTDIR)/save ]; then \
		mkdir -p $(OUTPUTDIR)/save; \
	fi
	@if [ ! -e $(OUTPUTDIR)/stat ]; then \
		mkdir -p $(OUTPUTDIR)/stat; \
	fi

datadel:
	$(RM) -r $(OUTPUTDIR)/log/*
	$(RM) -r $(OUTPUTDIR)/save/*
	$(RM) -r $(OUTPUTDIR)/stat/*

-include $(DEPS)

.PHONY : help all clean output datadel

