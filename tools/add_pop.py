"""Gabungkan populasi per sel (kolom pop_2026) ke grid.geojson sebagai properti 'pop' (jiwa, dibulatkan).
Pakai: python add_pop.py <penduduk_per_sel_riau_2026.csv> <data/grid.geojson>   (menimpa grid.geojson)"""
import json, sys, pandas as pd
csv, grid = sys.argv[1:3]
d = pd.read_csv(csv).set_index('grid_id')['pop_2026']
gj = json.load(open(grid))
miss = 0
for f in gj['features']:
    g = f['properties']['g']
    if g in d.index: f['properties']['pop'] = int(round(float(d[g])))
    else: f['properties']['pop'] = None; miss += 1
json.dump(gj, open(grid, 'w'), separators=(',', ':'))
print('sel', len(gj['features']), 'tanpa populasi', miss, 'total', sum(f['properties']['pop'] or 0 for f in gj['features']))
