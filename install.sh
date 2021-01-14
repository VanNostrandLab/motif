#!/usr/bin/env bash

set -e
unset PYTHONPATH

venv="$( cd .; pwd -P )"/venv
echo "Set up virtual environment for homer ...";

echo "    Installing Python (3.8) ...";
conda create --prefix="${venv}" --yes --quiet python=3.8 >/dev/null
echo "    Successful installed Python (3.8).";

echo "    Installing perl ...";
conda install --channel=bioconda --prefix="${venv}" --yes --quiet perl
echo "    Successful installed Perl.";

echo "    Installing homer ...";
conda install --channel=bioconda --prefix="${venv}" --yes --quiet homer >/dev/null
echo "    Successful installed homer."

echo "Successfully set up virtual environment for homer in ${venv}."