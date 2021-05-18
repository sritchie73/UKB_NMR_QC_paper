#!/bin/bash

outdir=data/raw/ukbiobank/extracted
indir=data/raw/ukbiobank/decoded

mkdir -p $outdir

./src/ukbtools/ukbconv $indir/ukb46936.enc_ukb csv \
  -e"$indir/encoding.ukb" \
  -o"$outdir/nightingale" \
  -i"$indir/fields.ukb"
mv fields.ukb $outdir

# Other formats for testing ukbnmr package

./src/ukbtools/ukbconv $indir/ukb46936.enc_ukb txt \
  -e"$indir/encoding.ukb" \
  -o"$outdir/nightingale" \
  -i"$indir/fields.ukb"

./src/ukbtools/ukbconv $indir/ukb46936.enc_ukb r \
  -e"$indir/encoding.ukb" \
  -o"$outdir/nightingale" \
  -i"$indir/fields.ukb"

./src/ukbtools/ukbconv $indir/ukb46936.enc_ukb docs \
  -e"$indir/encoding.ukb" \
  -o"$outdir/nightingale" \
  -i"$indir/fields.ukb"
rm fields.ukb

