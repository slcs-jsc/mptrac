for f in data/data.*/atm_2022_06_02_0*_00_00.tab; do
  ref="data.ref/${f#data/}"
  result=$(diff <(sort "$f") <(sort "$ref"))
  if [ -n "$result" ]; then
    echo "REAL DIFF: $f"
    echo "$result"
  fi
done
