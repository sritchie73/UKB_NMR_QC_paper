#!/bin/bash

outdir=data/raw/ukbiobank/extracted
indir=data/raw/ukbiobank/decoded

mkdir -p $outdir

./src/ukbtools/ukbconv $indir/ukb45386.enc_ukb csv \
  -e"$indir/encoding.ukb" \
  -o"$outdir/nightingale" \
  -i"$indir/fields.ukb"

