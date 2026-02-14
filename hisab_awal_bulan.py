import math
import sys
import os
import json
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
    jd_local = jd + tz_offset/24.0
    jdn_local = math.floor(jd_local + 0.5)
    f_local = (jd_local + 0.5 - jdn_local) * 24.0
    hours = int(f_local)
    minutes = int((f_local - hours) * 60)
    
    return f"{jdn_to_historical_str(jdn_local, False)} {hours:02d}:{minutes:02d}"

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
    
    def get_diff(jd):
        dt = get_delta_t(jd)
        # Sun position
        sun = VSOP87D_Sun.get_sun_position("VSOP87D.ear", jd + dt/86400.0)
        # Moon position (Geocentric Ecliptic)
        moon = elp.calculate(jd + dt/86400.0)
        
        diff = (moon['lon_deg'] - sun['lon_deg']) % 360
        if diff > 180: diff -= 360
        return diff

    # OPTIMASI 1: Wide enough search range (±2 days)
    current_jd = jd_start - 2.0
    prev_diff = get_diff(current_jd)
    
    target_jd = None
    # Scan with 2-hour steps over 4 days = 48 iterations
    for _ in range(48):
        current_jd += 2/24.0  # 2 hour steps
        diff = get_diff(current_jd)
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
                f_new = get_diff(x_new)
                
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
    def get_alt(jd):
        dt_sec = get_delta_t(jd)
        moon_res = elp.calculate_topocentric(jd, lat, lon, alt_obs_m=alt_m, delta_t_sec=dt_sec)
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
        alt = get_alt(curr)
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
                f_new = get_alt(x_new)
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

def draw_ascii_hilal(data):
    name = data['month']
    alt = data['alt']
    moon_az = data['moon_az']
    sun_az = data['sun_az']
    elong = data['elong']
    age = data['age']
    status = data['conclusion']
    greg_date = data['greg_date']
    sunset_str = data['sunset_time']
    moonset_str = data['moonset_time']
    sun_dist = data['sun_dist']
    moon_dist = data['moon_dist']
    tz_label = data['tz_label']
    est_date = data['est_date']
    next_m_name = data['month']

    print(f"\n{Color.BOLD}{Color.MAGENTA}--- VISUALISASI HILAL: {next_m_name.upper()} {greg_date} ---{Color.RESET}")
    # Colorize status
    status_col = Color.GREEN if status == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in status else Color.RED)
    print(f"Status: {status_col}{status}{Color.RESET}")
    print(f"Alt Hilal: {Color.GREEN}{alt:5.2f}\u00b0{Color.RESET} | Az Moon: {moon_az:6.2f}\u00b0 | Az Sun: {sun_az:6.2f}\u00b0 | Elong: {Color.YELLOW}{elong:5.2f}\u00b0{Color.RESET}")
    print(f"Ghurub   : {sunset_str} {tz_label} | Moonset: {moonset_str} {tz_label} | Umur: {age:5.2f} jam")
    print(f"Jarak Matahari: {sun_dist:8.6f} AU | Jarak Bulan: {moon_dist:7.0f} km")
    print(f"\n{Color.BOLD}KESIMPULAN: 1 {next_m_name} jatuh pada {Color.CYAN}{est_date}{Color.RESET}")
    
    if alt < -1:
        print(f"\n{Color.RED}[ Hilal jauh di bawah ufuk - Tidak terlihat ]{Color.RESET}")
        return

    # Dimensions
    rows = 22
    cols = 31
    rel_az_start = -15
    
    grid = [[" " for _ in range(cols)] for _ in range(rows)]
    
    # 1. Horizon line
    horizon_row = rows - 2
    for c in range(cols):
        grid[horizon_row][c] = f"{Color.GRAY}_{Color.RESET}"

    # 2. Mark Sun
    sun_col = 15 # center
    if 0 <= sun_col < cols:
        grid[horizon_row][sun_col] = f"{Color.YELLOW}S{Color.RESET}" # Sunset Point

    # 3. Mark Moon
    rel_az = moon_az - sun_az
    if rel_az > 180: rel_az -= 360
    if rel_az < -180: rel_az += 360
    
    moon_col = int(round(sun_col + rel_az))
    moon_row = horizon_row - int(round(alt * 2)) # 0.5 degrees per row
    
    is_visible_on_grid = False
    moon_char = f"{Color.WHITE}{Color.BOLD}({Color.RESET}"
    if 0 <= moon_row < horizon_row and 0 <= moon_col < cols:
        grid[moon_row][moon_col] = moon_char # Simple crescent
        is_visible_on_grid = True
    elif moon_row == horizon_row and 0 <= moon_col < cols:
        grid[moon_row][moon_col] = f"{Color.WHITE}C{Color.RESET}" # On horizon
        is_visible_on_grid = True

    # 4. Print Grid
    print(f"\n   {Color.GRAY}Alt (\u00b0){Color.RESET}")
    for r in range(rows - 1):
        if r % 4 == 0:
            alt_val = (horizon_row - r) / 2.0
            label = f"{Color.GRAY}{alt_val:4.1f} |{Color.RESET}"
        else:
            label = f"{Color.GRAY}     |{Color.RESET}"
        
        row_str = "".join(grid[r])
        print(f"{label}{row_str}")
    
    # X-axis label
    bottom_line = f"      {Color.GRAY}\u2534" + "\u2500" * (cols - 1) + f"{Color.RESET}"
    print(bottom_line)
    
    # Azimuth labels
    az_labels = f"      {Color.GRAY}"
    for i in range(rel_az_start, 16, 5):
        az_labels += f"{sun_az + i:3.0f}  "
    az_labels += Color.RESET
    print(az_labels)
    print(f"      {Color.GRAY}      (Azimut \u00b0){Color.RESET}")

    if not is_visible_on_grid:
        if alt < 0:
            print(f"\n{Color.RED}[ Note: Hilal di bawah horizon ]{Color.RESET}")
        else:
            print(f"\n{Color.YELLOW}[ Note: Hilal di luar jangkauan tampilan grid azimut ]{Color.RESET}")

    if status == "IMKAN RUKYAT":
        print(f"\n{Color.GREEN}{Color.BOLD}>>> Hilal kemungkinan besar terlihat (MABIMS terpenuhi) <<<{Color.RESET}")
    else:
        print(f"\n{Color.RED}>>> Hilal sulit/tidak terlihat (MABIMS tidak terpenuhi) <<<{Color.RESET}")

def run_hisab(hijri_year, lat, lon, alt_m, location_name, province_name=""):
    tz_offset, tz_label = get_political_tz(province_name, lon)
    
    print(f"\n{Color.BLUE}{Color.BOLD}===== HISAB AWAL BULAN {hijri_year} H ({location_name.upper()}) ====={Color.RESET}")
    print(f"Kriteria MABIMS: {Color.GREEN}Alt >= 3\u00b0{Color.RESET}, {Color.YELLOW}Elong >= 6.4\u00b0{Color.RESET} | Lat: {lat:.4f}, Lon: {lon:.4f}, Alt: {alt_m}m")
    print(f"Zona Waktu: {Color.CYAN}{tz_label}{Color.RESET} (UTC+{tz_offset})")
    print(f"{Color.GRAY}" + "-" * 150 + f"{Color.RESET}")
    header = f"{'No':<2} | {'Bulan':<15} | {f'Ijtima ({tz_label})':<16} | {f'Ghurub ({tz_label})':<16} | {'Umur':<10} | {'Alt':<7} | {'Elong':<7} | {'Status':<15} | {'Kesimpulan'}"
    print(f"{Color.BOLD}{header}{Color.RESET}")
    print(f"{Color.GRAY}" + "-" * 150 + f"{Color.RESET}")
    
    monthly_results = []

    def calculate_night(jd_target_sunset, jd_conj, ijtima_str, label_prefix="", month_idx=0):
        # Convert JD to historical date string
        jdn_sunset = math.floor(jd_target_sunset + 0.5)
        
        # Time part in local TZ
        jd_local = jd_target_sunset + tz_offset/24.0
        f_local = (jd_local + 0.5 - math.floor(jd_local + 0.5)) * 24.0
        h_local, m_local = int(f_local), int((f_local - int(f_local)) * 60)
        sunset_str = f"{jdn_to_historical_str(jdn_sunset, False)} {h_local:02d}:{m_local:02d}"

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

        # Colorize result
        c_alt = Color.GREEN if moon_res['alt_apparent'] >= 3.0 else (Color.YELLOW if moon_res['alt_apparent'] >= 0 else Color.RED)
        c_elong = Color.GREEN if elong >= 6.4 else Color.YELLOW
        c_conclusion = Color.GREEN if conclusion == "IMKAN RUKYAT" else (Color.YELLOW if "ISTIKMAL" in conclusion else Color.RED)

        # Number or H+1
        lbl_num = f"{month_idx+1:2}" if "(H+1)" not in label_prefix else "  "
        
        row = (f"{lbl_num} | {label_prefix:<15} | {ijtima_str:<16} | {sunset_str:<16} | "
               f"{age_hours:<8.2f} j | {c_alt}{moon_res['alt_apparent']:<6.2f}\u00b0{Color.RESET} | "
               f"{c_elong}{elong:<6.2f}\u00b0{Color.RESET} | {status:<15} | {c_conclusion}{Color.BOLD}{conclusion}{Color.RESET}")
        print(row)
        
        if "(H+1)" not in label_prefix:
            # Estimate 1st day of month
            days_to_add = 1 if conclusion == "IMKAN RUKYAT" else 2
            est_jdn = jdn_sunset + days_to_add
            est_date_str = jdn_to_historical_str(est_jdn)

            jd_moonset = find_moonset(jd_target_sunset, lat, lon, alt_m)
            if jd_moonset:
                jd_ms_local = jd_moonset + tz_offset/24.0
                f_ms = (jd_ms_local + 0.5 - math.floor(jd_ms_local + 0.5)) * 24.0
                ms_str = f"{int(f_ms):02d}:{int((f_ms-int(f_ms))*60):02d}"
            else:
                ms_str = "--:--"

            monthly_results.append({
                'month': HIJRI_MONTH_NAMES[month_idx],
                'greg_date': jdn_to_historical_str(jdn_sunset),
                'alt': moon_res['alt_apparent'],
                'moon_az': moon_res['az_deg'],
                'sun_az': topo_sun['azimuth'],
                'sunset_time': f"{h_local:02d}:{m_local:02d}",
                'moonset_time': ms_str,
                'tz_label': tz_label,
                'est_date': est_date_str,
                'sun_dist': sun_pos['range'],
                'moon_dist': moon_geo['dist_km'],
                'elong': elong,
                'age': age_hours,
                'conclusion': conclusion
            })
            
        return trigger_h1

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
        is_late = calculate_night(jd_sunset_1, jd_ijtima, ijtima_str, HIJRI_MONTH_NAMES[m-1], m-1)
        
        # If Ijtima was after sunset, check the next night too
        if is_late:
            jd_sunset_2 = jd_sunset_1 + 1.0
            calculate_night(jd_sunset_2, jd_ijtima, ijtima_str, f"  (H+1)", m-1)
            
        print() # Line break for readability
    
    return monthly_results

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
    
    # Prepare the formatted strings (with index)
    # items should be strings like " 1. Aceh"
    
    # Calculate max width for each column
    max_len = max(len(s) for s in items) + 2
    
    for r in range(rows):
        line = ""
        for c in range(cols):
            idx = r + c * rows
            if idx < len(items):
                line += items[idx].ljust(max_len)
        print(line)

def select_location(provinces, cities):
    print("\n--- PILIH PROVINSI ---")
    prov_list = [f"{i+1:2}. {p['name']}" for i, p in enumerate(provinces)]
    print_grid(prov_list, cols=3)
    
    choice = input("\nMasukkan nomor provinsi (atau 'keluar'): ").strip().lower()
    if choice == 'keluar': return "EXIT"
    
    try:
        idx = int(choice) - 1
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
    
    choice = input("\nMasukkan nomor lokasi (atau 'kembali' / 'keluar'): ").strip().lower()
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
            continue # Try again if invalid or none picked
            
        year_str = input(f"\nLokasi: {loc['type']} {loc['name']}, {loc['province']}\nMasukkan Tahun Hijriyah (e.g. 1447) atau 'keluar': ").strip().lower()
        if year_str == 'keluar':
            break
            
        try:
            hijri_year = int(year_str)
            results = run_hisab(hijri_year, loc['lat'], loc['lon'], loc['elevation_m'], f"{loc['type']} {loc['name']}", loc['province'])
            
            while True:
                viz_choice = input("\nMasukkan nomor bulan untuk visualisasi ASCII (1-12) atau Enter untuk kembali: ").strip()
                if not viz_choice: break
                if viz_choice.lower() == 'keluar': return
                
                try:
                    m_idx = int(viz_choice) - 1
                    if 0 <= m_idx < len(results):
                        draw_ascii_hilal(results[m_idx])
                    else:
                        print("Nomor bulan tidak valid (1-12).")
                except ValueError:
                    print("Masukkan angka 1-12.")
                    
        except ValueError:
            print("Tahun tidak valid.")
            
        input("\nTekan Enter untuk kembali ke awal...")

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
            print("\nProgram dihentikan.")
            sys.exit(0)
