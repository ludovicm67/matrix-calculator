# We define variables for directories
DIROBJ  = obj
DIRINC  = include
DIRSRC  = src
DIRBIN  = bin

# Compiler
CC      = gcc -g
CFLAGS  = -W -Wall
LDFLAGS = -lm

# Dependencies, objects, ...
DEPS    = $(wildcard include/*.h)
SOURCES = $(wildcard src/*.c)
OBJETS  = $(SOURCES:src/%.c=obj/%.o)
EXEC    = main

.SUFFIXES:

# We create targets
$(DIRBIN)/$(EXEC) : $(OBJETS) | $(DIRBIN)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(DIROBJ)/%.o : $(DIRSRC)/%.c $(DEPS) | $(DIROBJ)
	$(CC) $(CFLAGS) -c -I$(DIRINC) -o $@ $<

$(DIROBJ) $(DIRBIN):
	@mkdir -p $@

# Targets to call manually
.PHONY: archive
archive:
	tar -f archive.tar.gz -cvz $(DIRSRC)/*.c $(DIRINC)/*.h Makefile

.PHONY: clean
clean:
	rm -rf $(DIROBJ) $(DIRBIN)
