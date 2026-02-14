import math
import sys
import os
import json
import textwrap
from datetime import datetime, timedelta

# Append path to ensure modules are found
sys.path.append(os.getcwd())

from matahari import VSOP87D_Sun, EarthRotation
from bulan import ELP2000_Moon

# --- ANSI COLORS ---
class Color:
    CYAN = "\033[96m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BLUE = "\033[94m"
    MAGENTA = "\033[95m"
    WHITE = "\033[97m"
    GRAY = "\033[90m"
    BOLD = "\033[1m"
    RESET = "\033[0m"

    @staticmethod
    def rgb(r, g, b):
        return f"\033[38;2;{r};{g};{b}m"

    @staticmethod
    def bg_rgb(r, g, b):
        return f"\033[48;2;{r};{g};{b}m"

# --- KUWAITI ALGORITHM (Ported from JS) ---

HIJRI_EPOCH_JDN = 1948440
CYCLE_DAYS = 10631
CYCLE_YEARS = 30
LEAP_YEARS_KUWAITI = [2, 5, 7, 10, 13, 15, 18, 21, 24, 26, 29]
HIJRI_MONTH_NAMES = [
    "Muharram", "Safar", "Rabiul Awal", "Rabiul Akhir",
    "Jumadil Awal", "Jumadil Akhir", "Rajab", "Sya'ban",
    "Ramadhan", "Syawal", "Dzulkaidah", "Dzulhijjah"
]

def is_hijri_leap(year):
    cycle_year = ((year - 1) % 30) + 1
    return cycle_year in LEAP_YEARS_KUWAITI

def get_hijri_year_length(year):
    return 355 if is_hijri_leap(year) else 354

def get_hijri_month_lengths(year):
    lengths = [30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29]
    if is_hijri_leap(year):
        lengths[11] = 30
    return lengths

def hijri_to_jdn(day, month, year):
    ii = year - 1
    cycles = ii // 30
    year_in_cycle = ii % 30
    
    days = cycles * 10631
    for k in range(1, year_in_cycle + 1):
        days += 354
        if k in LEAP_YEARS_KUWAITI:
            days += 1
            
    month_lengths = get_hijri_month_lengths(year)
    for m in range(month - 1):
        days += month_lengths[m]
        
    days += day
    return days + HIJRI_EPOCH_JDN - 1

def jdn_to_gregorian(jdn):
    l = jdn + 68569
    n = math.floor((4 * l) / 146097)
    l1 = l - math.floor((146097 * n + 3) / 4)
    i = math.floor((4000 * (l1 + 1)) / 1461001)
    l2 = l1 - math.floor((1461 * i) / 4) + 31
    j = math.floor((80 * l2) / 2447)
    day = l2 - math.floor((2447 * j) / 80)
    l3 = math.floor(j / 11)
    month = j + 2 - 12 * l3
    year = 100 * (n - 49) + i + l3
    return day, month, year

def jdn_to_julian(jdn):
    j = jdn + 1402
    k = (j - 1) // 1461
    l = j - 1461 * k
    m = (l - 1) // 31
    n = (l - 1) // 153
    p = (n * 2 + 1) // 5
    day = l - 153 * n + p * 153 // 2 # Simplified but usually JDN conversion is standard
    # Using Meeus formula for Julian:
    b = jdn + 1524
    c = math.floor((b - 122.1) / 365.25)
    d = math.floor(365.25 * c)
    e = math.floor((b - d) / 30.6001)
    day = b - d - math.floor(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    return day, month, year

def jdn_to_historical_str(jdn, include_calendar=True):
    """Converts JDN to a string using Julian before 1582-10-15 and Gregorian after."""
    if jdn < 2299161: # Before 15 Oct 1582 (Gregorian jump)
        d, m, y = jdn_to_julian(jdn)
        suffix = " (Jul)" if include_calendar else ""
    else:
        d, m, y = jdn_to_gregorian(jdn)
        suffix = " (Greg)" if include_calendar else ""
    return f"{d:02d}-{m:02d}-{y}{suffix}"

def jdn_to_jd(jdn):
    return jdn - 0.5

def get_political_tz(province_name, lon):
    """
    Returns (offset_hours, label) based on Indonesian political time zones.
    If outside Indonesia, fallback to longitude-based offset.
    """
    p = province_name.upper()
    
    # WIB: Sumatra, Java, West/Central Kalimantan
    wib_provinces = [
        "ACEH", "SUMATERA UTARA", "SUMATERA BARAT", "RIAU", "JAMBI", 
        "SUMATERA SELATAN", "BENGKULU", "LAMPUNG", "KEPULAUAN BANGKA BELITUNG", 
        "KEPULAUAN RIAU", "DKI JAKARTA", "JAWA BARAT", "JAWA TENGAH", 
        "DI YOGYAKARTA", "JAWA TIMUR", "BANTEN", "KALIMANTAN BARAT", 
        "KALIMANTAN TENGAH"
    ]
    
    # WITA: Sulawesi, Bali, NTB, NTT, South/East/North Kalimantan
    wita_provinces = [
        "BALI", "NUSA TENGGARA BARAT", "NUSA TENGGARA TIMUR", 
        "KALIMANTAN SELATAN", "KALIMANTAN TIMUR", "KALIMANTAN UTARA", 
        "SULAWESI UTARA", "SULAWESI TENGAH", "SULAWESI SELATAN", 
        "SULAWESI TENGGARA", "GORONTALO", "SULAWESI BARAT"
    ]
    
    # WIT: Maluku, Papua
    wit_provinces = [
        "MALUKU", "MALUKU UTARA", "PAPUA BARAT", "PAPUA", 
        "PAPUA SELATAN", "PAPUA TENGAH", "PAPUA PEGUNUNGAN", "PAPUA BARAT DAYA"
    ]

    if p in wib_provinces: return 7, "WIB"
    if p in wita_provinces: return 8, "WITA"
    if p in wit_provinces: return 9, "WIT"
    
    # Fallback for international or unknown
    offset = round(lon / 15.0)
    label = f"UTC{'+' if offset >= 0 else ''}{offset}"
    return offset, label

def jd_to_historical_time_str(jd, tz_offset, tz_label):
    """Converts JD to a string with date and time using specific TZ."""
    if jd is None: return "--"
    
    jd_local = jd + tz_offset/24.0
    jdn_local = math.floor(jd_local + 0.5)
    f_local = (jd_local + 0.5 - jdn_local) * 24.0
    
    # Precise hour/min/sec calculation
    total_seconds = int(round(f_local * 3600))
    hours = (total_seconds // 3600) % 24
    minutes = (total_seconds // 60) % 60
    seconds = total_seconds % 60
    
    # Handle JN rollover if hours reached 24 from rounding
    if total_seconds >= 86400:
        jdn_local += 1
        hours = 0
    
    return f"{jdn_to_historical_str(jdn_local, False)} {hours:02d}:{minutes:02d}:{seconds:02d}"

# --- ASTRONOMICAL ENGINE ---

elp = ELP2000_Moon(data_dir=os.getcwd())
elp.load_all_series()

def get_delta_t(jd):
    """
    Kalkulasi Delta-T (\u0394T) berbasis polinomial Espenak & Meeus (NASA).
    Akurasi tinggi untuk rentang sejarah luas.
    """
    y = 2000 + (jd - 2451545.0) / 365.2425
    
    if y < -500:
        u = (y - 1820) / 100
        return -20 + 32 * u**2
    elif y < 500:
        u = y / 100
        return 10583.6 - 1014.41*u + 33.78311*u**2 - 5.952053*u**3 - 0.1798452*u**4 + 0.022174192*u**5 + 0.0090316521*u**6
    elif y < 1600:
        u = (y - 1000) / 100
        return 1574.2 - 556.01*u + 71.23472*u**2 + 0.319781*u**3 - 0.8503463*u**4 - 0.005050998*u**5 + 0.0083572073*u**6
    elif y < 1700:
        u = (y - 1600) / 100
        return 120 - 98.3*u - 153.2*u**2 + 163 * u**3
    elif y < 1800:
        u = (y - 1700) / 100
        return 8.83 + 160.3*u - 59.2*u**2 + 133.3*u**3 - 14.74*u**4
    elif y < 1860:
        u = (y - 1800) / 100
        return 13.72 - 33.28*u + 273.05*u**2 - 27.8181*u**3 + 1.2287*u**4 - 0.0161*u**5
    elif y < 1900:
        u = (y - 1860) / 100
        return 7.62 + 57.37*u - 251.75*u**2 + 16.8066*u**3 - 0.44736*u**4 + 0.00505*u**5
    elif y < 1920:
        u = (y - 1900) / 100
        return -2.79 + 26.83*u + 11.23*u**2 - 2.012*u**3
    elif y < 1941:
        u = (y - 1920) / 100
        return 21.20 + 84.49*u - 205.3*u**2 + 12.89*u**3
    elif y < 1961:
        u = (y - 1950) / 100
        return 29.07 + 40.7*u + 13.0*u**2 - 22.99*u**3
    elif y < 1986:
        u = (y - 1975) / 100
        return 45.45 + 106.7*u - 4.1*u**2 - 3.7*u**3
    elif y < 2005:
        u = (y - 2000) / 100
        return 63.86 + 33.488 * u - 269.8 * u**2 + 674.3 * u**3
    elif y < 2050:
        u = y - 2000
        return 62.92 + 0.32217 * u + 0.005589 * u**2
    else:
        u = (y - 1820) / 100
        return -20 + 32 * u**2

def get_conjunction(jd_start):
    """
    Finds the exact moment of conjunction (Ijtima) around jd_start (JD UT1).
    
    OPTIMASI: Menggunakan Secant Method + Smart Initial Guess
    - Lebih cepat converge dari bisection
    - Gunakan estimasi synodic month untuk initial guess
    """
    delta_t_ref = get_delta_t(jd_start)
    
    def get_diff(jd, apparent=True):
        dt = get_delta_t(jd)
        # Sun position (Apparent)
        sun = VSOP87D_Sun.get_sun_position("VSOP87D.ear", jd + dt/86400.0, apparent=apparent)
        # Moon position (Geocentric Ecliptic Apparent)
        moon = elp.calculate(jd + dt/86400.0, apparent=apparent)
        
        diff = (moon['lon_deg'] - sun['lon_deg']) % 360
        if diff > 180: diff -= 360
        return diff

    # OPTIMASI 1: Wide enough search range (±2 days)
    current_jd = jd_start - 2.0
    prev_diff = get_diff(current_jd)
    
    target_jd = None
    # Scan with 2-hour steps over 4 days = 48 iterations
    # OPTIMASI: Use apparent=False for fast scanning
    for _ in range(48):
        current_jd += 2/24.0  # 2 hour steps
        diff = get_diff(current_jd, apparent=False)
        if (prev_diff < 0 and diff >= 0) or (prev_diff > 0 and diff <= 0):
            # Zero crossing found! Gunakan Secant Method untuk converge lebih cepat
            x0 = current_jd - 2/24.0
            x1 = current_jd
            f0 = prev_diff
            f1 = diff
            
            # Secant Method - typically converges in 5-8 iterations
            for _ in range(10):
                if abs(f1 - f0) < 1e-10:
                    break
                # Secant: x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
                x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
                f_new = get_diff(x_new, apparent=True) # Full precision for refinement
                
                x0, f0 = x1, f1
                x1, f1 = x_new, f_new
                
                # Jika sudah sangat dekat, stop
                if abs(f1) < 0.001:  # 0.001 degree ~ 3.6 arcsec
                    break
            
            target_jd = x1
            break
        prev_diff = diff
    
    return target_jd

def find_moonset(jd_sunset, lat, lon, alt_m=10):
    """Finds Moonset JD around the time of Sunset.
    
    OPTIMASI: Reduce iterations + use Secant Method
    """
    def get_alt(jd, apparent=True):
        dt_sec = get_delta_t(jd)
        moon_res = elp.calculate_topocentric(jd, lat, lon, alt_obs_m=alt_m, delta_t_sec=dt_sec, apparent=apparent)
        # Moon semidiameter ~0.25 + refraction ~0.583 total
        return moon_res['alt_apparent'] + 0.583 
    
    # OPTIMASI: Reduce search window from ±4h to ±3h
    low_limit = jd_sunset - 3/24.0
    high_limit = jd_sunset + 3/24.0
    
    curr = low_limit
    step = 1/24.0  # 1 hour steps (lebih efisien)
    prev_alt = get_alt(curr)
    
    while curr < high_limit:
        curr += step
        # OPTIMASI: Use apparent=False for scanning
        alt = get_alt(curr, apparent=False)
        if prev_alt > 0 and alt <= 0:
            # Moonset found! Use Secant Method
            x0 = curr - step
            x1 = curr
            f0 = prev_alt
            f1 = alt
            
            # Secant Method
            for _ in range(8):
                if abs(f1 - f0) < 1e-10:
                    break
                x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
                f_new = get_alt(x_new, apparent=True) # Full precision for refinement
                x0, f0 = x1, f1
                x1, f1 = x_new, f_new
                
                if abs(f1) < 0.1:  # 0.1 degree threshold
                    break
            return x1
        prev_alt = alt
    return None

def find_sunset(jd_date, lat, lon, alt_m=10, tz_offset=7):
    """Finds Sunset JD (UT1) on the Gregorian day of jd_date (in local time).
    
    OPTIMASI: Reduce iterations + use Secant Method
    """
    # 1. Get the local date
    jd_local_mid = jd_date + tz_offset/24.0
    jdn_local = math.floor(jd_local_mid + 0.5)
    
    # Approx JD for 12:00 local time (UT)
    jd_mid_ut = jdn_local - tz_offset/24.0
    
    def get_alt(jd):
        dt_sec = get_delta_t(jd)
        sun_pos = VSOP87D_Sun.get_sun_position("VSOP87D.ear", jd + dt_sec/86400.0)
        ra, dec = EarthRotation.ecliptic_to_equatorial(sun_pos['lon'], sun_pos['lat'], jd)
        topo = EarthRotation.to_topocentric(ra, dec, jd, lat, lon, alt_obs_m=alt_m)
        return topo['altitude_apparent'] + 0.8333 

    # OPTIMASI: Reduce search range
    low = jd_mid_ut + 4/24.0
    high = jd_mid_ut + 8/24.0
    
    # Get initial values for Secant Method
    f_low = get_alt(low)
    f_high = get_alt(high)
    
    # Check if sunset exists in range
    if f_low <= 0 or f_high >= 0:
        # Use bisection as fallback (for extreme latitudes)
        for _ in range(15):
            mid = (low + high) / 2
            if get_alt(mid) > 0:
                low = mid
            else:
                high = mid
        return (low + high) / 2
    
    # Secant Method
    for _ in range(12):
        if abs(f_high - f_low) < 1e-10:
            break
        x_new = high - f_high * (high - low) / (f_high - f_low)
        f_new = get_alt(x_new)
        
        low, f_low = high, f_high
        high, f_high = x_new, f_new
        
        if abs(f_high) < 0.01:  # 0.01 degree threshold
            break
    
    return high

def show_logo():
    os.system('cls' if os.name == 'nt' else 'clear')
    print(f"""{Color.CYAN}{Color.BOLD}
  ██╗  ██╗ ████████╗ ███████╗
  ██║  ██║ ╚══██╔══╝ ██╔════╝
  ███████║    ██║    █████╗  
  ██╔══██║    ██║    ██╔══╝  
  ██║  ██║    ██║    ██║     
  ╚═╝  ╚═╝    ╚═╝    ╚═╝     {Color.WHITE}
  Hilal Tracker Falakiyah
{Color.GRAY}  |  Kreator: Tri Setyo Pamungkas
  |  Engine : VSOP87 & ELP2000-82
  |  Website: hisab.tawheed.my.id{Color.RESET}
    """)

def draw_ascii_hilal(data_list):
    # data_list can be a single dict or a list of two dicts (Day 29 and H+1)
    if isinstance(data_list, dict):
        data_list = [data_list]
    
    is_dual = len(data_list) > 1
    
    # Header Info (based on the first day, but we'll show both in plots)
    data = data_list[0]
    next_m_name = data['month']
    greg_date = data['greg_date']
    status = data['conclusion']
    est_date = data['est_date']

    print(f"\n{Color.BOLD}{Color.MAGENTA}--- VISUALISASI HILAL (SIDE-BY-SIDE) : {next_m_name.upper()} {greg_date} ---{Color.RESET}")
    
    # 1. Prepare Grids for each day
    grids = []
    meta_infos = []
    
    rows = 16 # Reduced from 20 to make it less "tall"
    cols = 61
    horizon_row = rows - 3
    sun_col = cols // 2
    
    # Calculate Dynamic Scaling
    max_alt = max([d['alt'] for d in data_list])
    # Standard scale: 1.5 rows per degree
    # Available sky rows: horizon_row (0 to horizon_row-1) = rows-3
    available_sky_rows = horizon_row
    
    # We want max_alt to fit at row 1 (leaving row 0 for margin)
    # Target row index = 1
    # moon_row = horizon_row - (alt * scale)
    # 1 = horizon_row - (max_alt * scale)
    # scale = (horizon_row - 1) / max_alt
    
    if max_alt > (available_sky_rows - 1) / 1.5:
        # If max_alt is too high for standard scale, reduce scale
        alt_scale = (available_sky_rows - 1) / max_alt
    else:
        alt_scale = 1.5 # Standard
    
    for d_idx, d in enumerate(data_list):
        alt = d['alt']
        moon_az = d['moon_az']
        sun_az = d['sun_az']
        
        # Sky Gradient
        sky_colors = []
        for r in range(horizon_row):
            p = r / (horizon_row - 1) if horizon_row > 1 else 1.0
            p_color = p**1.5
            red = int(20 * (1-p) + 255 * p_color)
            green = int(20 * (1-p) + 120 * p_color)
            blue = int(80 * (1-p) + 20 * p_color)
            sky_colors.append(Color.bg_rgb(red, green, blue))

        grid = [[sky_colors[r] + " " + Color.RESET for c in range(cols)] for r in range(horizon_row)]
        
        # Stars
        import random
        random.seed(42 + d_idx) # slightly different stars per day
        for r in range(int(horizon_row * 0.6)):
            for _ in range(2):
                c = random.randint(0, cols - 1)
                if r < horizon_row * 0.4:
                    star_color = Color.WHITE if random.random() > 0.3 else Color.GRAY
                    grid[r][c] = sky_colors[r] + star_color + ('.' if random.random() > 0.5 else '*') + Color.RESET

        # Horizon
        horizon_str = Color.GREEN + "▔" + Color.RESET
        grid.append([horizon_str for _ in range(cols)])
        grid.append([Color.GREEN + "░" + Color.RESET for _ in range(cols)])
        grid.append([Color.GREEN + " " + Color.RESET for _ in range(cols)])

        # Sun
        if 0 <= sun_col < cols:
            grid[horizon_row][sun_col] = Color.YELLOW + "●" + Color.RESET

        # Moon
        rel_az = moon_az - sun_az
        if rel_az > 180: rel_az -= 360
        if rel_az < -180: rel_az += 360
        
        moon_col = int(round(sun_col + rel_az * 2))
        moon_row = horizon_row - int(round(alt * alt_scale))
        
        moon_char = "☽" if rel_az > 0 else "☾"
        if alt > -2: # show moon even if slightly below horizon but within grid
            if 0 <= moon_row < horizon_row and 1 <= moon_col < cols - 1:
                grid[moon_row][moon_col] = sky_colors[moon_row] + Color.WHITE + Color.BOLD + moon_char + Color.RESET
            elif moon_row == horizon_row and 1 <= moon_col < cols - 1:
                 grid[moon_row][moon_col] = Color.WHITE + Color.BOLD + moon_char + Color.RESET

        # Azimuth Labels & Markers
        labels_row = [" " for _ in range(cols)]
        marker_row = [" " for _ in range(cols)]
        min_az = sun_az - 15
        max_az = sun_az + 15
        start_mark = math.ceil(min_az / 5.0) * 5
        for az_val in range(int(start_mark), int(max_az) + 1, 5):
            c = int(round(sun_col + (az_val - sun_az) * 2))
            if 0 <= c <= cols - 3:
                val_str = f"{az_val % 360:3.0f}"
                for i, char in enumerate(val_str):
                    if c + i < cols: labels_row[c + i] = char
        
        if 0 <= sun_col < cols: marker_row[sun_col] = f"{Color.YELLOW}^{Color.RESET}"
        if 0 <= moon_col < cols: marker_row[moon_col] = f"{Color.WHITE}{Color.BOLD}^{Color.RESET}"

        grids.append({
            'main': grid,
            'labels': "".join(labels_row),
            'markers': "".join(marker_row),
            'date': d['greg_date'],
            'day_lbl': "Hari ke-29" if d_idx == 0 else "H+1 (Hari ke-30)",
            'alt': d['alt']
        })

    # 2. Print Side-by-Side
    spacing = "    "
    
    # Sub-headers
    sub_header = ""
    for g in grids:
        lbl = f"({g['day_lbl']}: {g['date']})"
        sub_header += f"{Color.BOLD}{lbl:^67}{Color.RESET}{spacing}"
    print(f"\n{sub_header}")
    
    # Detailed Info Rows
    # We need to construct lines that have info for grid 1, spacing, info for grid 2
    
    # Extract data for easier access
    infos = []
    for d in data_list:
        status_c = Color.GREEN if d['conclusion'] == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in d['conclusion'] else Color.RED)
        infos.append({
            'status_str': f"Status: {status_c}{d['conclusion']}{Color.RESET}",
            'pos_str': f"Alt: {Color.GREEN}{d['alt']:5.2f}\\u00b0{Color.RESET} | Az M: {d['moon_az']:6.2f}\\u00b0 | Az S: {d['sun_az']:6.2f}\\u00b0",
            'elong_str': f"Elong: {Color.YELLOW}{d['elong']:5.2f}\\u00b0{Color.RESET} | Umur: {d['age']:5.2f} jam",
            'time_str': f"Ghurub: {d['sunset_time']} | Moonset: {d['moonset_time']} {d['tz_label']}",
            'dist_str': f"S.Dist: {d['sun_dist']:6.4f} AU | M.Dist: {d['moon_dist']:6.0f} km"
        })
        
    for key in ['status_str', 'pos_str', 'elong_str', 'time_str', 'dist_str']:
        line_str = ""
        for i, info in enumerate(infos):
            # We need to pad visual length to 67
            text = info[key]
            # Simple approach: just assume it fits or use fixed spacing if 1 grid
            if i < len(infos) - 1:
                # Pad based on visible length? Too complex for quick script.
                # Let's just use tab or fixed large padding?
                # The grid is 67 chars wide.
                # Let's try to align by appending spaces.
                # Hacky:
                visible_len = len(text.replace(Color.GREEN, "").replace(Color.YELLOW, "").replace(Color.RED, "").replace(Color.RESET, "").replace(Color.BOLD, "").replace(Color.CYAN, ""))
                padding = 67 - visible_len
                line_str += text + " " * padding + spacing
            else:
                line_str += text
        print(line_str)

    print("-" * (67 * len(grids) + 4 * (len(grids)-1)))
    
    print(f" Alt | {' '*59} |{spacing} Alt | {' '*59} |")
    for r in range(rows):
        line = ""
        for g_idx, g in enumerate(grids):
            if r < horizon_row and r % 4 == 0:
                alt_val = (horizon_row - r) / alt_scale
                prefix = f"{alt_val:4.1f} | "
            elif r == horizon_row:
                prefix = f" 0.0 | "
            else:
                prefix = f"     | "
            
            line += prefix + "".join(g['main'][r]) + spacing
        print(line)

    # X-axis
    x_line = ""
    for g in grids:
        x_line += f" Alt | {Color.GRAY}\u2514" + "\u2500" * (cols - 1) + f"{Color.RESET}{spacing}"
    print(x_line)

    # Azimuth Labels
    az_line = ""
    for g in grids:
        az_line += f"  Az | {g['labels']}{spacing}"
    print(az_line)

    # Markers
    mk_line = ""
    for g in grids:
        mk_line += f"     | {g['markers']}{spacing}"
    print(mk_line)
    
    print(f"     | {Color.YELLOW}S=Matahari{Color.RESET}  {Color.WHITE}{Color.BOLD}M=Bulan{Color.RESET}             {spacing}     | {Color.YELLOW}S=Matahari{Color.RESET}  {Color.WHITE}{Color.BOLD}M=Bulan{Color.RESET}")

    # Final Summary
    print(f"\n{Color.BOLD}KESIMPULAN: 1 {next_m_name} jatuh pada {Color.CYAN}{est_date}{Color.RESET}")
    if status == "IMKAN RUKYAT":
        print(f"{Color.GREEN}{Color.BOLD}>>> HASIL: Hilal TERLIHAT pada Hari ke-29 (MABIMS terpenuhi) <<<{Color.RESET}")
    else:
        print(f"{Color.RED}{Color.BOLD}>>> HASIL: Hilal TIDAK TERLIHAT pada Hari ke-29 -> ISTIKMAL <<<{Color.RESET}")
        print(f"{Color.RED}{Color.BOLD}>>> 1 {next_m_name} jatuh pada hari berikutnya (H+1) <<<{Color.RESET}")

def calculate_night_single(jd_target_sunset, jd_conj, ijtima_str, lat, lon, alt_m, tz_offset, tz_label, label_prefix="", month_idx=0):
    # Convert JD to historical date string
    jdn_sunset = math.floor(jd_target_sunset + 0.5)
    
    # Time part in local TZ
    jd_local = jd_target_sunset + tz_offset/24.0
    f_local = (jd_local + 0.5 - math.floor(jd_local + 0.5)) * 24.0
    h_local = int(f_local)
    m_local = int((f_local - h_local) * 60)
    s_local = int(round(((f_local - h_local) * 60 - m_local) * 60))
    if s_local == 60:
        m_local += 1
        s_local = 0
    if m_local == 60:
        h_local += 1
        m_local = 0
        
    sunset_str = f"{jdn_to_historical_str(jdn_sunset, False)} {h_local:02d}:{m_local:02d}:{s_local:02d}"

    age_hours = (jd_target_sunset - jd_conj) * 24.0
    ijtima_after_sunset = (jd_conj > jd_target_sunset)
    
    # Calculation logic
    sun_pos = VSOP87D_Sun.get_sun_position("VSOP87D.ear", jd_target_sunset + get_delta_t(jd_target_sunset)/86400.0)
    ra_sun, dec_sun = EarthRotation.ecliptic_to_equatorial(sun_pos['lon'], sun_pos['lat'], jd_target_sunset)
    topo_sun = EarthRotation.to_topocentric(ra_sun, dec_sun, jd_target_sunset, lat, lon, alt_obs_m=alt_m)
    moon_geo = elp.calculate(jd_target_sunset + get_delta_t(jd_target_sunset)/86400.0)
    moon_res = elp.calculate_topocentric(jd_target_sunset, lat, lon, alt_obs_m=alt_m, delta_t_sec=get_delta_t(jd_target_sunset))

    cos_elong = (math.sin(moon_geo['lat_rad'])*math.sin(sun_pos['lat']) + 
                 math.cos(moon_geo['lat_rad'])*math.cos(sun_pos['lat'])*math.cos(moon_geo['lon_rad'] - sun_pos['lon']))
    elong = math.degrees(math.acos(max(-1.0, min(1.0, cos_elong))))

    if ijtima_after_sunset:
        status, conclusion, trigger_h1 = "BELUM IJTIMA", "MUSTAHIL", True
    elif moon_res['alt_apparent'] < 0:
        status, conclusion, trigger_h1 = "BELUM UFUK", "MUSTAHIL", True
    else:
        status, trigger_h1 = "SUDAH UFUK", False
        mabims_met = (moon_res['alt_apparent'] >= 3.0) and (elong >= 6.4)
        conclusion = "IMKAN RUKYAT" if mabims_met else "ISTIKMAL (30)"

    res_data = {
        'month_idx': month_idx,
        'label_prefix': label_prefix,
        'ijtima_str': ijtima_str,
        'sunset_str': sunset_str,
        'h_local': h_local, 'm_local': m_local,
        'age': age_hours,
        'alt': moon_res['alt_apparent'],
        'elong': elong,
        'status': status,
        'conclusion': conclusion,
        'moon_az': moon_res['az_deg'],
        'sun_az': topo_sun['azimuth'],
        'sun_dist': sun_pos['range'],
        'moon_dist': moon_geo['dist_km'],
        'lat': lat, 'lon': lon, 'tz_label': tz_label, 'tz_offset': tz_offset,
        'jdn_sunset': jdn_sunset,
        'jd_target_sunset': jd_target_sunset
    }
    
    # Calculate estimates
    if "(H+1)" not in label_prefix:
        days_to_add = 1 if conclusion == "IMKAN RUKYAT" else 2
        est_jdn = jdn_sunset + days_to_add
        res_data['est_date'] = jdn_to_historical_str(est_jdn)
        res_data['est_jdn'] = est_jdn
    else:
        est_jdn = jdn_sunset + 1
        res_data['est_date'] = jdn_to_historical_str(est_jdn)
        res_data['est_jdn'] = est_jdn
        
    # Calculate moonset
    jd_moonset = find_moonset(jd_target_sunset, lat, lon, alt_m)
    if jd_moonset:
        jd_ms_local = jd_moonset + tz_offset/24.0
        f_ms = (jd_ms_local + 0.5 - math.floor(jd_ms_local + 0.5)) * 24.0
        res_data['moonset_time'] = f"{int(f_ms):02d}:{int((f_ms-int(f_ms))*60):02d}"
    else:
        res_data['moonset_time'] = "--:--"
    
    res_data['trigger_h1'] = trigger_h1
    return res_data

def run_hisab(hijri_year, lat, lon, alt_m, location_name, province_name=""):
    tz_offset, tz_label = get_political_tz(province_name, lon)
    
    print(f"\n{Color.BLUE}{Color.BOLD}===== HISAB AWAL BULAN {hijri_year} H ({location_name.upper()}) ====={Color.RESET}")
    print(f"Kriteria MABIMS: {Color.GREEN}Alt >= 3\u00b0{Color.RESET}, {Color.YELLOW}Elong >= 6.4\u00b0{Color.RESET} | Lat: {lat:.4f}, Lon: {lon:.4f}, Alt: {alt_m}m")
    print(f"Zona Waktu: {Color.CYAN}{tz_label}{Color.RESET} (UTC+{tz_offset})")
    print(f"{Color.GRAY}" + "-" * 155 + f"{Color.RESET}")
    header = f"{'No':<2} | {'Bulan':<15} | {f'Ijtima ({tz_label})':<19} | {f'Ghurub ({tz_label})':<19} | {'Umur':<10} | {'Alt':<7} | {'Elong':<7} | {'Status':<15} | {'Kesimpulan'}"
    print(f"{Color.BOLD}{header}{Color.RESET}")
    print(f"{Color.GRAY}" + "-" * 155 + f"{Color.RESET}")
    
    monthly_results = []
    
    for m in range(1, 13):
        # Anchor on Day 29 of previous month
        prev_m = m - 1
        prev_y = hijri_year
        if prev_m == 0:
            prev_m = 12
            prev_y -= 1
        
        jdn_29 = hijri_to_jdn(29, prev_m, prev_y)
        jd_anchor = jdn_to_jd(jdn_29)
        
        jd_ijtima = get_conjunction(jd_anchor)
        if not jd_ijtima: continue
        
        ijtima_str = jd_to_historical_time_str(jd_ijtima, tz_offset, tz_label)

        # Night 1: The day of Ijtima
        jd_sunset_1 = find_sunset(jd_ijtima, lat, lon, alt_m, tz_offset)
        
        # Calculate Night 1
        d1 = calculate_night_single(jd_sunset_1, jd_ijtima, ijtima_str, lat, lon, alt_m, tz_offset, tz_label, HIJRI_MONTH_NAMES[m-1], m-1)
        
        # Print Row 1
        c_alt = Color.GREEN if d1['alt'] >= 3.0 else (Color.YELLOW if d1['alt'] >= 0 else Color.RED)
        c_elong = Color.GREEN if d1['elong'] >= 6.4 else Color.YELLOW
        c_conclusion = Color.GREEN if d1['conclusion'] == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in d1['conclusion'] else Color.RED)
        
        lbl_num = f"{m:2}"
        row = (f"{lbl_num} | {d1['label_prefix']:<15} | {d1['ijtima_str']:<19} | {d1['sunset_str']:<19} | "
               f"{d1['age']:<8.2f} j | {c_alt}{d1['alt']:<6.2f}\u00b0{Color.RESET} | "
               f"{c_elong}{d1['elong']:<6.2f}\u00b0{Color.RESET} | {d1['status']:<15} | {c_conclusion}{Color.BOLD}{d1['conclusion']}{Color.RESET}")
        print(row)
        
        # Prepare result entry
        res_entry = {
            'month': HIJRI_MONTH_NAMES[m-1],
            'greg_date': jdn_to_historical_str(d1['jdn_sunset']),
            'alt': d1['alt'],
            'moon_az': d1['moon_az'],
            'sun_az': d1['sun_az'],
            'sunset_time': f"{d1['h_local']:02d}:{d1['m_local']:02d}",
            'moonset_time': d1['moonset_time'],
            'tz_label': d1['tz_label'],
            'est_date': d1['est_date'],
            'sun_dist': d1['sun_dist'],
            'moon_dist': d1['moon_dist'],
            'elong': d1['elong'],
            'age': d1['age'],
            'conclusion': d1['conclusion']
        }
        monthly_results.append(res_entry)

        # If Ijtima was after sunset, check the next night too
        if d1['trigger_h1']:
            jd_sunset_2 = jd_sunset_1 + 1.0
            d2 = calculate_night_single(jd_sunset_2, jd_ijtima, ijtima_str, lat, lon, alt_m, tz_offset, tz_label, f"  (H+1)", m-1)
            
            # Print Row 2 (H+1)
            c_alt = Color.GREEN if d2['alt'] >= 3.0 else (Color.YELLOW if d2['alt'] >= 0 else Color.RED)
            c_elong = Color.GREEN if d2['elong'] >= 6.4 else Color.YELLOW
            c_conclusion = Color.GREEN if d2['conclusion'] == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in d2['conclusion'] else Color.RED)
            
            row2 = (f"   | {d2['label_prefix']:<15} | {d2['ijtima_str']:<19} | {d2['sunset_str']:<19} | "
                   f"{d2['age']:<8.2f} j | {c_alt}{d2['alt']:<6.2f}\u00b0{Color.RESET} | "
                   f"{c_elong}{d2['elong']:<6.2f}\u00b0{Color.RESET} | {d2['status']:<15} | {c_conclusion}{Color.BOLD}{d2['conclusion']}{Color.RESET}")
            print(row2)
            
            # Store H+1 data
            data_h1 = {
                'month': HIJRI_MONTH_NAMES[m-1],
                'greg_date': jdn_to_historical_str(d2['jdn_sunset']),
                'alt': d2['alt'],
                'moon_az': d2['moon_az'],
                'sun_az': d2['sun_az'],
                'sunset_time': f"{d2['h_local']:02d}:{d2['m_local']:02d}",
                'moonset_time': d2['moonset_time'],
                'tz_label': d2['tz_label'],
                'est_date': d2['est_date'],
                'sun_dist': d2['sun_dist'],
                'moon_dist': d2['moon_dist'],
                'elong': d2['elong'],
                'age': d2['age'],
                'conclusion': d2['conclusion']
            }
            res_entry['h1_data'] = data_h1
            
        print() # Line break for readability
    
    return monthly_results

def run_hisab_all_indonesia(hijri_year, month_int, cities):
    m_idx = month_int - 1
    m_name = HIJRI_MONTH_NAMES[m_idx]
    
    print(f"\n{Color.BLUE}{Color.BOLD}===== HISAB SELURUH KOTA/KAB INDONESIA: {m_name.upper()} {hijri_year} H ====={Color.RESET}")
    # Determine next month name for header
    next_m_idx = (m_idx + 1) % 12
    next_m_name = HIJRI_MONTH_NAMES[next_m_idx]
    col_date_header = f"1 {next_m_name}"

    print(f"{Color.GRAY}" + "-" * 175 + f"{Color.RESET}") # Extended width
    header = f"{'No':<4} | {'Kota/Kabupaten':<30} | {'Provinsi':<20} | {'Ijtima (Local)':<14} | {'Ghurub':<8} | {'Alt':<7} | {'Elong':<7} | {'Status':<15} | {'Kesimpulan':<15} | {col_date_header}"
    print(f"{Color.BOLD}{header}{Color.RESET}")
    print(f"{Color.GRAY}" + "-" * 175 + f"{Color.RESET}")
    
    # Calculate Ijtima (Global, but need anchor)
    prev_m = month_int - 1
    prev_y = hijri_year
    if prev_m == 0:
        prev_m = 12
        prev_y -= 1
        
    jdn_29 = hijri_to_jdn(29, prev_m, prev_y)
    jd_anchor = jdn_to_jd(jdn_29)
    jd_ijtima = get_conjunction(jd_anchor)
    
    if not jd_ijtima:
        print("Ijtima tidak ditemukan.")
        return

    days_id = ["Ahad", "Senin", "Selasa", "Rabu", "Kamis", "Jumat", "Sabtu"]

    # Loop cities
    for i, c in enumerate(cities):
        city_name = f"{c['type']} {c['name']}"
        prov_name = c['province']
        lat = c['lat']
        lon = c['lon']
        alt_m = c['elevation_m'] if 'elevation_m' in c else 10
        
        tz_offset, tz_label = get_political_tz(prov_name, lon)
        ijtima_str = jd_to_historical_time_str(jd_ijtima, tz_offset, tz_label)
        ijtima_time = ijtima_str.split(' ')[1] if ' ' in ijtima_str else ijtima_str
        
        jd_sunset = find_sunset(jd_ijtima, lat, lon, alt_m, tz_offset)
        d = calculate_night_single(jd_sunset, jd_ijtima, ijtima_str, lat, lon, alt_m, tz_offset, tz_label, m_name, m_idx)
        
        # Colorize
        c_alt = Color.GREEN if d['alt'] >= 3.0 else (Color.YELLOW if d['alt'] >= 0 else Color.RED)
        c_elong = Color.GREEN if d['elong'] >= 6.4 else Color.YELLOW
        c_conclusion = Color.GREEN if d['conclusion'] == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in d['conclusion'] else Color.RED)
        
        # Calculate full date string
        est_jdn = d.get('est_jdn', 0)
        day_idx = int((est_jdn + 1.5) % 7)
        full_date = f"{days_id[day_idx]}, {d['est_date']}"

        # Wrap Text
        city_lines = textwrap.wrap(city_name, width=30)
        prov_lines = textwrap.wrap(prov_name, width=20)
        max_lines = max(len(city_lines), len(prov_lines))
        
        for j in range(max_lines):
            c_txt = city_lines[j] if j < len(city_lines) else ""
            p_txt = prov_lines[j] if j < len(prov_lines) else ""
            
            if j == 0:
                row = (f"{i+1:<4} | {c_txt:<30} | {p_txt:<20} | {ijtima_time:<14} | {d['sunset_str'][-8:]:<8} | "
                       f"{c_alt}{d['alt']:<6.2f}\u00b0{Color.RESET} | "
                       f"{c_elong}{d['elong']:<6.2f}\u00b0{Color.RESET} | {d['status']:<15} | {c_conclusion}{Color.BOLD}{d['conclusion']:<15}{Color.RESET} | {full_date}")
            else:
                row = (f"{'':<4} | {c_txt:<30} | {p_txt:<20} | {'':<14} | {'':<8} | "
                       f"{'':<7} | "
                       f"{'':<7} | {'':<15} | {'':<15} | ")
            print(row)
    print(f"{Color.GRAY}" + "-" * 175 + f"{Color.RESET}")


def load_data():
    with open('provinces.json', 'r') as f:
        provinces = json.load(f)
    with open('cities_indonesia_processed.json', 'r') as f:
        cities = json.load(f)
    return provinces, cities

def print_grid(items, cols=3):
    """Prints a list of strings in a clean grid format."""
    if not items: return
    rows = math.ceil(len(items) / cols)
    
    # Max width
    max_len = max(len(s) for s in items) + 2
    
    for r in range(rows):
        line = ""
        for c in range(cols):
            idx = r + c * rows
            if idx < len(items):
                line += items[idx].ljust(max_len)
        print(line)

def select_location(provinces, cities):
    print(f"\n--- PILIH PROVINSI ---")
    print(" 0. HISAB SELURUH KOTA/KAB INDONESIA") # Option 0
    prov_list = [f"{i+1:2}. {p['name']}" for i, p in enumerate(provinces)]
    print_grid(prov_list, cols=3)
    
    choice = input(f"\nMasukkan nomor provinsi (0-34) (atau 'keluar'): ").strip().lower()
    if choice == 'keluar': return "EXIT"
    
    try:
        idx = int(choice)
        if idx == 0: return "ALL"
        
        idx -= 1
        if not (0 <= idx < len(provinces)): raise ValueError
        selected_prov = provinces[idx]['name']
    except (ValueError, IndexError):
        print("Pilihan tidak valid.")
        return select_location(provinces, cities)

    # Filter cities by province
    filtered_locations = [c for c in cities if c['province'].upper() == selected_prov.upper()]
    
    if not filtered_locations:
        print(f"Tidak ada data kota/kabupaten untuk provinsi {selected_prov}.")
        return select_location(provinces, cities)

    print(f"\n--- KOTA/KABUPATEN DI {selected_prov} ---")
    loc_list = [f"{i+1:3}. {c['type'][:3]}. {c['name']}" for i, c in enumerate(filtered_locations)]
    print_grid(loc_list, cols=3 if len(loc_list) < 30 else 4)
    
    choice = input(f"\nMasukkan nomor lokasi (atau 'kembali' / 'keluar'): ").strip().lower()
    if choice == 'keluar': return "EXIT"
    if choice == 'kembali': return select_location(provinces, cities)

    try:
        idx = int(choice) - 1
        if not (0 <= idx < len(filtered_locations)): raise ValueError
        return filtered_locations[idx]
    except (ValueError, IndexError):
        print("Pilihan tidak valid.")
        return select_location(provinces, cities)

def main():
    provinces, cities = load_data()
    
    while True:
        os.system('cls' if os.name == 'nt' else 'clear')
        show_logo()
        
        loc = select_location(provinces, cities)
        if loc == "EXIT" or loc is None:
            if loc == "EXIT": break
            continue 
            
        year_str = input(f"\nMasukkan Tahun Hijriyah (e.g. 1447) atau 'keluar': ").strip().lower()
        if year_str == 'keluar':
            break
            
        try:
            hijri_year = int(year_str)
            if loc == "ALL":
                # Special flow for All Indonesia
                print(f"\n--- PILIH BULAN ---")
                month_list = [f"{i+1:2}. {n}" for i, n in enumerate(HIJRI_MONTH_NAMES)]
                print_grid(month_list, cols=3)
                m_choice = input(f"\nMasukkan nomor bulan (1-12): ").strip()
                try:
                    m_int = int(m_choice)
                    if 1 <= m_int <= 12:
                         run_hisab_all_indonesia(hijri_year, m_int, cities)
                    else:
                        print("Bulan tidak valid.")
                except ValueError:
                    print("Input bulan salah.")
                input(f"\nTekan Enter untuk kembali ke awal...")
            else:
                results = run_hisab(hijri_year, loc['lat'], loc['lon'], loc['elevation_m'], f"{loc['type']} {loc['name']}", loc['province'])
                
                while True:
                    viz_choice = input(f"\nMasukkan nomor bulan untuk visualisasi ASCII (1-12) atau Enter untuk kembali: ").strip()
                    if not viz_choice: break
                    if viz_choice.lower() == 'keluar': return
                    
                    try:
                        m_idx = int(viz_choice) - 1
                        if 0 <= m_idx < len(results):
                            res = results[m_idx]
                            if 'h1_data' in res:
                                draw_ascii_hilal([res, res['h1_data']])
                            else:
                                draw_ascii_hilal(res)
                        else:
                            print("Nomor bulan tidak valid (1-12).")
                    except ValueError:
                        print("Masukkan angka 1-12.")
                    
        except ValueError:
            print("Tahun tidak valid.")
            input("Tekan Enter untuk kembali...")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # CLI usage (retaining functionality for backward compatibility/automation)
        try:
            h_year = int(sys.argv[1])
            # Default to Jakarta if run via CLI without interactive pick
            run_hisab(h_year, -6.2088, 106.8456, 12, "Jakarta")
        except ValueError:
            print("Usage: python hisab_awal_bulan.py [HIJRI_YEAR]")
    else:
        try:
            main()
        except KeyboardInterrupt:
            print(f"\nProgram dihentikan.")
            sys.exit(0)
