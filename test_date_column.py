
import sys
import os
sys.path.append(os.getcwd())
# Redirect stdout
sys.stdout = open('test_output_date.txt', 'w', encoding='utf-8')

from hisab_awal_bulan import load_data, run_hisab_all_indonesia

print("Testing Date Column in All Indonesia Hisab")
provinces, cities = load_data()

# Pick specific cities
test_cities = [c for c in cities if c['name'] in ["Banda Aceh", "Jakarta", "Merauke"]]
if test_cities:
    run_hisab_all_indonesia(1447, 9, test_cities)
else:
    print("Test cities not found.")
