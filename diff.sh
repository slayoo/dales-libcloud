#/bin/bash
diff -ruN \
  --exclude=build \
  --exclude=modlibcloud.f90 \
  dales.orig dales \
  > dales.diff
