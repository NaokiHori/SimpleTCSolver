CC       := cc
CFLAGS   := -std=c99 -Wall -Wextra -O3
DEPEND   := -MMD
LIBS     := -lm
INCLUDES := -Iinclude
SRCSDIR  := src
OBJSDIR  := obj
SRCS     := $(foreach dir, $(shell find $(SRCSDIR) -type d), $(wildcard $(dir)/*.c))
OBJS     := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.o)))
DEPS     := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.d)))
TARGET   := a.out

help:
	@echo "all     : create \"$(TARGET)\""
	@echo "clean   : remove \"$(TARGET)\" and object files \"$(OBJSDIR)/*.o\""
	@echo "output  : make directory to store NPY files"
	@echo "datadel : remove NPY files"
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
	@if [ ! -e output ]; then \
		mkdir output; \
	fi

datadel:
	$(RM) output/*.npy

-include $(DEPS)

.PHONY : help all clean output datadel

