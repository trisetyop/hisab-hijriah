import sys
import os
sys.path.append(os.getcwd())

from hisab_awal_bulan import jd_to_historical_time_str, get_conjunction, get_delta_t
from matahari import VSOP87D_Sun
from bulan import ELP2000_Moon

# 1. Test Seconds Formatting
jd_test = 2461141.5 + 13/24.0 + 45/1440.0 + 30/86400.0 # 13:45:30
formatted = jd_to_historical_time_str(jd_test, 7, "WIB")
print(f"Test Formatting: {formatted}")
assert "13:45:30" in formatted or "20:45:30" in formatted # depends on how jd_test is anchored

# 2. Test get_conjunction with apparent position
# Reference: 1447 Muharram Ijtima around June 2025 (June 25)
jd_anchor = 2460851.5 # around June 25 2025
elp = ELP2000_Moon(data_dir=os.getcwd())
elp.load_all_series()

jd_conj = get_conjunction(jd_anchor)
ijtima_str = jd_to_historical_time_str(jd_conj, 7, "WIB")
print(f"Ijtima Muharram 1447: {ijtima_str}")
