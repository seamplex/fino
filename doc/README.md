This directory contains the Fino documentation.
Everything is created out of the source files which are a mixture of M4 and Markdown.

# Sources

Run `make.sh`

## Manpage

Edit `fino.1.md` 


## Manual

Edit `fino.md`.
Pay attention to `reference-manual.md`.


## Help screen

Edit `help.md`.

**TODO:** read it from the commented sources.

# Syntax highlighting files

## Kate

Run

```
./syntax-kate.sh > $HOME/.kde/share/apps/katepart/syntax/fino.xml
./syntax-kate.sh > $HOME/.local/share/katepart5/syntax/fino.xml
```

## TeX

Generate `keywords.tex` with

```
$ ./syntax-tex.sh > keywords.tex
```

and then include `syntax.tex` (which includes `keywords.tex`) in your preamble

## Geany (unmaintained)

Run

```
./syntax-geany.sh > $HOME/.config/geany/filedefs/filetypes.Fino.conf
```

Remember to add the line

```
Fino=*.fin;
```

o `filetype_extensions.conf` (tools -> configuration files -> filetype_extensions.conf)

