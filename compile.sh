#!/usr/bin/bash

latexmk -cd ./informe/document -pdf -latexoption=-file-line-error -latexoption=-interaction=nonstopmode && cp informe/document.pdf TP3-G01.pdf
