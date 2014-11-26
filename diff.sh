#/bin/bash
diff -ruN \
  --exclude=build \
  dales.orig dales \
  > dales.diff
