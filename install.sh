#!/usr/bin/env bash

set -e
unset PYTHONPATH

venv="$( cd .; pwd -P )"/venv
echo "Set up virtual environment for motif pipeline ...";

if [ ! -f "${venv}/bin/python" ]
  echo "    Installing Python (3.8) ...";
  conda create --prefix="${venv}" --yes --quiet python=3.8 >/dev/null
  echo "    Successful installed Python (3.8)."
fi

if [ ! -f "${venv}/bin/perl" ]
  echo "    Installing perl ...";
  conda install --channel=bioconda --prefix="${venv}" --yes --quiet perl
  echo "    Successful installed Perl."
fi

if [ ! -f "${venv}/bin/python" ]
  echo "    Installing homer ...";
  conda install --channel=bioconda --prefix="${venv}" --yes --quiet homer >/dev/null
  echo "    Successful installed homer."
fi

echo "    Finalizing installation ...";
cp "${venv}/source/annotate_peaks_bedformat_with_proxdistal_lncRNA.pl" "${venv}/bin/"
chmod +x "${venv}/bin/annotate_peaks_bedformat_with_proxdistal_lncRNA.pl"

cp "${venv}/source/pull_seqs_for_regions_and_run_homer_clip_analysis_version.pl" "${venv}/bin/"
chmod +x "${venv}/bin/pull_seqs_for_regions_and_run_homer_clip_analysis_version.pl"

cp "${venv}/source/motif.py" "${venv}/bin/"

read -r -d '' script << EOF || true
#!/usr/bin/env bash

unset PYTHONPATH
export TMPDIR=/storage/vannostrand/tmp
export TEMP=/storage/vannostrand/tmp
export TMP=/storage/vannostrand/tmp

export PATH="${venv}/bin:$PATH

python ${venv}/bin/motif.py \$@
EOF

echo "${script}" > "${clip}/motif"
chmod +x "${clip}/motif"
echo "    Successfully finalized installation."

echo "Successfully set up virtual environment for motif pipeline in ${venv}."