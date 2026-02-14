import urllib.request
import os
import math
from functools import lru_cache

file_url = "https://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87D.ear"
file_name = "VSOP87D.ear"

if not os.path.exists(file_name):
    print(f"Downloading {file_name}...")
    urllib.request.urlretrieve(file_url, file_name)
    print("Download selesai.")

class VSOP87D_Sun:
    """
    Parser khusus VSOP87D untuk menghitung posisi Matahari Geosentrik
    berdasarkan data Helioselektrik Bumi (VSOP87D.ear).
    
    OPTIMASI: Menggunakan caching untuk menyimpan data file yang sudah di-parse.
    """

    DPI = 2.0 * math.pi
    T2000 = 2451545.0
    A1000 = 365250.0
    
    # ========== CACHE UNTUK DATA FILE VSOP87D ==========
    _vsop_cache = {}  # Cache untuk data file yang sudah di-parse

    @classmethod
    def get_sun_position(cls, filepath, jd, prec=0.0):
        """
        Menghitung posisi Matahari (Geosentrik) dari file Bumi (VSOP87D.ear).

        Returns:
            dict: {
                'lon': Bujur Matahari (radian),
                'lat': Lintang Matahari (radian),
                'range': Jarak (AU),
                'lon_deg': Bujur (derajat)
            }
        """
        # 1. Hitung posisi Helioselektrik Bumi
        # ivers=4 (VSOP87D), ibody=3 (Earth)
        r_earth = cls._compute_vsop87d(filepath, jd, 4, "EARTH", prec)

        # 2. Konversi ke Geosentrik Matahari
        # L_sun = L_earth + PI
        lon_sun = (r_earth[0] + math.pi) % cls.DPI

        # B_sun = -B_earth
        lat_sun = -r_earth[1]

        # R tetap sama
        dist = r_earth[2]

        return {
            'lon': lon_sun,
            'lat': lat_sun,
            'range': dist,
            'lon_deg': math.degrees(lon_sun)
        }

    @classmethod
    def _parse_vsop87_file(cls, filepath, ivers, body_name):
        """
        Parsing file VSOP87D dan simpan ke cache.
        Setiap file hanya di-parse SEKALI saja.
        """
        cache_key = (filepath, ivers, body_name)
        
        # Jika sudah ada di cache, return yang sudah di-cache
        if cache_key in cls._vsop_cache:
            return cls._vsop_cache[cache_key]
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File data tidak ditemukan: {filepath}")
        
        # Parse file dan simpan ke list
        terms = []  # List of (ic, it, a, b, c)
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        idx = 0
        while idx < len(lines):
            line = lines[idx]
            if len(line) < 67: break
            
            # Header format
            iv = int(line[17:18])
            bo = line[22:29].strip()
            ic = int(line[41:42])
            it = int(line[59:60])
            n_terms = int(line[60:67])
            idx += 1
            
            # Skip jika bukan Bumi atau versi D
            if iv != ivers or bo != body_name:
                idx += n_terms
                continue
            
            # Parse semua terms
            for _ in range(n_terms):
                term_line = lines[idx]
                idx += 1
                
                a = float(term_line[79:97])
                b = float(term_line[97:111])
                c = float(term_line[111:131])
                
                terms.append((ic, it, a, b, c))
        
        # Simpan ke cache
        cls._vsop_cache[cache_key] = terms
        print(f"[CACHE] VSOP87D data loaded: {len(terms)} terms for {body_name}")
        
        return terms
    
    @classmethod
    def _compute_vsop87d(cls, filepath, tdj, ivers, body_name, prec):
        """
        Internal parser untuk membaca file VSOP87D.
        OPTIMASI: Menggunakan data dari cache, tidak perlu parsing ulang.
        """
        # Ambil data dari cache (otomatis di-parse jika belum ada)
        terms = cls._parse_vsop87_file(filepath, ivers, body_name)
        
        # Hitung T (millennia dari J2000.0)
        t_powers = [1.0] * 7
        t_powers[1] = (tdj - cls.T2000) / cls.A1000
        for i in range(2, 7):
            t_powers[i] = t_powers[1] * t_powers[i-1]
        
        r = [0.0] * 6
        
        # Hitung menggunakan data dari cache
        for ic, it, a, b, c in terms:
            # Evaluasi A * cos(B + C*T) * T^k
            u = b + c * t_powers[1]
            r[ic - 1] += a * math.cos(u) * t_powers[it]
        
        return r

def julian_now():
    """Mendapatkan Julian Date saat ini (UT) - Sederhana."""
    from datetime import datetime, timezone
    now = datetime.now(timezone.utc)
    # Rumus praktis sederhana untuk JD
    Y, M, D = now.year, now.month, now.day
    if M <= 2:
        Y -= 1
        M += 12
    A = math.floor(Y / 100)
    B = 2 - A + math.floor(A / 4)
    JD = math.floor(365.25 * (Y + 4716)) + math.floor(30.6001 * (M + 1)) + D + B - 1524.5
    # Tambah fraksi jam
    day_fraction = (now.hour + now.minute/60 + now.second/3600) / 24
    return JD + day_fraction

# --- CONTOH EKSEKUSI ---
if __name__ == "__main__":
    # Pastikan file VSOP87D.ear ada di folder yang sama
    data_file = "VSOP87D.ear"

    jd_input = julian_now()
    print(f"Menghitung untuk Julian Date: {jd_input:.6f}")

    try:
        sun = VSOP87D_Sun.get_sun_position(data_file, jd_input)

        print("\n--- HASIL POSISI MATAHARI (GEOSENTRIK) ---")
        print(f"Bujur Ekliptika (λ) : {sun['lon_deg']:.6f}°")
        print(f"Lintang Ekliptika (β): {math.degrees(sun['lat']):.10f}°")
        print(f"Jarak (R)           : {sun['range']:.8f} AU")

        # Info tambahan untuk Ilmu Falak
        # Menghitung perkiraan kasar letak buruj
        buruj = ["Hamal", "Tsaur", "Jauza", "Sarthon", "Asad", "Sunbulah",
                 "Mizan", "Aqrob", "Qous", "Jadyu", "Dalwu", "Hut"]
        idx_buruj = int(sun['lon_deg'] / 30) % 12
        print(f"Posisi Buruj        : {buruj[idx_buruj]}")

    except FileNotFoundError:
        print(f"Error: File '{data_file}' tidak ditemukan.")
        print("Silakan unduh file VSOP87D.ear dari sumber resmi Bureau des Longitudes.")



import math

class EarthRotation:
    """
    Modul rotasi bumi, sidereal time, dan konversi koordinat
    Ekliptika -> Ekuatorial -> Toposentrik (Azimut/Alt).
    """

    @staticmethod
    def get_obliquity(jd):
        """Mean Obliquity of the Ecliptic (Laskar / Meeus)"""
        T = (jd - 2451545.0) / 36525.0
        eps0 = (23.43929111 -
                (46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3) / 3600.0)
        return math.radians(eps0)

    @staticmethod
    def get_dip(elevation_m):
        """
        Menghitung Kerendahan Ufuk (Dip of Horizon).
        Rumus: 0.0293 * sqrt(h_meter) derajat.
        """
        if elevation_m <= 0: return 0.0
        return 0.0293 * math.sqrt(elevation_m)

    @staticmethod
    def get_refraction(alt_deg, pressure=1010, temp_c=10):
        """
        Menghitung Refraksi Atmosfer (Meeus / Bennett).
        Input: Altitude geometris (derajat).
        Output: Koreksi refraksi (derajat).
        Apparent Alt = Geometric Alt + Refraction.
        """
        # Refraksi tidak berlaku di bawah horizon jauh (-5 derajat?)
        # Rumus Bennett tidak stabil jika h dekat -90.
        h = max(alt_deg, -5.0)
        
        # R dalam arcmin
        R_arcmin = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
        
        # Koreksi Suhu (C) dan Tekanan (mbar)
        R_val = R_arcmin * (pressure / 1010.0) * (283.0 / (273.0 + temp_c))
        
        return R_val / 60.0 # derajat

    @staticmethod
    def get_nutation(jd):
        """
        Menghitung Nutasi (IAU 1980 Reduced Accuracy / Meeus Ch. 22).
        Returns:
            (delta_psi, delta_epsilon) dalam derajat.
            delta_psi: Nutasi Longitude
            delta_eps: Nutasi Obliquity
        """
        T = (jd - 2451545.0) / 36525.0
        
        # Argumen Fundamental (derajat)
        # D: Elongasi Bulan
        D = (297.85036 + 445267.111480 * T - 0.0019142 * T**2 + T**3/189474) % 360
        # M: Anomali Matahari
        M = (357.52772 + 35999.050340 * T - 0.0001603 * T**2 - T**3/300000) % 360
        # M': Anomali Bulan
        Mp = (134.96298 + 477198.867398 * T + 0.0086972 * T**2 + T**3/56250) % 360
        # F: Argumen Lintang Bulan
        F = (93.27191 + 483202.017538 * T - 0.0036825 * T**2 + T**3/327270) % 360
        # Omega: Longitude Node Bulan
        Om = (125.04452 - 1934.136261 * T + 0.0020708 * T**2 + T**3/450000) % 360
        
        D_r = math.radians(D)
        M_r = math.radians(M)
        Mp_r = math.radians(Mp)
        F_r = math.radians(F)
        Om_r = math.radians(Om)
        
        # Suku-suku utama Nutasi (arcseconds) - Meeus Table 22.A
        # Term: (d, m, m', f, om) -> (coef_sin_psi, coef_cos_eps)
        # 1. -17.20 sin(Om)
        d_psi = -17.20 * math.sin(Om_r)
        d_eps = 9.20 * math.cos(Om_r)
        
        # 2. -1.32 sin(2L) -> L ~ F + Om ? No, L = Mean Longitude Sun
        # L_sun ~ 280.4665 + 36000.7698*T
        # But Meeus args are D, M, M', F, Omega.
        # Main terms from low accuracy list:
        
        # -1.32 sin(2Lsun) -> -1.32 sin(2(F+Om+...)) ??
        # Let's use the exact lines from Meeus Low Precision:
        # Longitude (d_psi):
        d_psi += -1.32 * math.sin(2*F_r - 2*D_r + 2*Om_r) # sin(2L_sun) ?? Check definition.
        # Actually standard series 2L is 2*(280.47 + 36000.77T).
        # Which corresponds to 2*(F - D + Omega)? No.
        # L = F + Omega + D? No.
        # Let's use standard args directly:
        
        # -1.32 sin(2L) : 2L = 2*MeanLonSun. MeanLonSun = 280.4665 + ...
        # Or derived: 2L = 2*F - 2*D + 2*Omega + 2*M? No.
        # Let's stick to using the 4 main terms usually cited (Meeus pg 144):
        # 1. Om: -17.20 sin(Om)
        # 2. 2L: -1.32 sin(2L)  (L = Mean Longitude Sun)
        # 3. 2L': -0.23 sin(2L') (L' = Mean Longitude Moon)
        # 4. 2Om: +0.21 sin(2Om)
        
        # Need L and Lp.
        L = (280.4665 + 36000.7698 * T) % 360
        L_r = math.radians(L)
        
        Lp = (218.3165 + 481267.8813 * T) % 360
        Lp_r = math.radians(Lp)
        
        d_psi += -1.32 * math.sin(2 * L_r)
        d_psi += -0.23 * math.sin(2 * Lp_r)
        d_psi += 0.21 * math.sin(2 * Om_r)
        
        # Obliquity (d_eps):
        d_eps += 0.57 * math.cos(2 * L_r)
        d_eps += 0.10 * math.cos(2 * Lp_r)
        d_eps += -0.09 * math.cos(2 * Om_r)
        
        return (d_psi / 3600.0), (d_eps / 3600.0) # degrees

    @staticmethod
    def get_aberration_low_prec(dist_au):
        """
        Menghitung Konstanta Aberasi (Approximation).
        Delta Lambda = -20.49" / R
        Delta Beta = 0
        """
        val_sec = -20.4898 / dist_au
        return val_sec / 3600.0

    @staticmethod
    def get_gmst(jd_ut):
        """Greenwich Mean Sidereal Time (derajat)"""
        T = (jd_ut - 2451545.0) / 36525.0
        gmst = (280.46061837 +
                360.98564736629 * (jd_ut - 2451545.0) +
                0.000387933 * T**2 -
                (T**3 / 38710000.0))
        return gmst % 360.0

    @classmethod
    def ecliptic_to_equatorial(cls, lon_rad, lat_rad, jd, nut_lon_deg=0.0, nut_obl_deg=0.0):
        """
        Konversi koordinat ekliptika (λ, β) -> ekuatorial (RA, Dec)
        Input dari VSOP87 (spherical).
        
        Jika nut_lon_deg & nut_obl_deg != 0, menghitung 'True Equator & Equinox of Date'.
        lon_rad += nut_lon (konversi ke rad dulu).
        Obliquity = Mean + nut_obl.
        """
        eps = cls.get_obliquity(jd) + math.radians(nut_obl_deg)
        
        # Apply Nutation in Longitude to the input Longitude
        lon_true = lon_rad + math.radians(nut_lon_deg)

        sin_ra = (math.sin(lon_true) * math.cos(eps) -
                  math.tan(lat_rad) * math.sin(eps))
        cos_ra = math.cos(lon_true)
        ra_rad = math.atan2(sin_ra, cos_ra)
        ra_rad = (ra_rad + 2 * math.pi) % (2 * math.pi)

        sin_dec = (math.sin(lat_rad) * math.cos(eps) +
                   math.cos(lat_rad) * math.sin(eps) * math.sin(lon_true))
        dec_rad = math.asin(sin_dec)

        return ra_rad, dec_rad

    @classmethod
    def to_topocentric(cls, ra_rad, dec_rad, jd_ut, lat_obs, lon_obs, alt_obs_m=0.0, pressure=1010, temp_c=10):
        """
        Konversi RA/Dec -> Azimut & Altitude (toposentrik geometrik)
        
        lat_obs, lon_obs dalam derajat
        Azimut: 0° = Utara, 90° = Timur
        Output termasuk 'altitude_apparent' (refraksi) dan 'altitude_observed' (+dip).
        """

        # Local Sidereal Time
        gmst = cls.get_gmst(jd_ut)
        lst_rad = math.radians((gmst + lon_obs) % 360.0)

        # Hour Angle
        h = lst_rad - ra_rad

        phi = math.radians(lat_obs)
        delta = dec_rad

        # Altitude
        sin_alt = (math.sin(phi) * math.sin(delta) +
                   math.cos(phi) * math.cos(delta) * math.cos(h))
        alt_rad = math.asin(sin_alt)
        alt_deg = math.degrees(alt_rad)

        # Azimuth (Meeus formula 13.5)
        # 0 = North, 90 = East, 180 = South, 270 = West
        # tan A = sin(H) / (cos(H)sin(phi) - tan(delta)cos(phi))
        # With atan2(x, y):
        # x = -sin(H) * cos(delta)
        # y = sin(delta)*cos(phi) - cos(delta)*cos(H)*sin(phi)
        
        sin_h = math.sin(h)
        cos_h = math.cos(h)
        sin_dec = math.sin(delta)
        cos_dec = math.cos(delta)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)
        
        x = -sin_h * cos_dec
        y = (sin_dec * cos_phi) - (cos_dec * cos_h * sin_phi)
        
        az_rad = math.atan2(x, y)
        
        # Atmospheric
        dip = cls.get_dip(alt_obs_m)
        refr = cls.get_refraction(alt_deg - dip, pressure, temp_c)
        
        alt_apparent = alt_deg + refr
        alt_observed = alt_apparent + dip # Relative to Sea Horizon if Dip > 0

        return {
            "azimuth": (math.degrees(az_rad) + 360.0) % 360.0,
            "altitude": alt_deg,
            "altitude_apparent": alt_apparent,
            "altitude_observed": alt_observed,
            "dip": dip,
            "refraction": refr
        }



# Note: Examples moved to __main__ block




import math
from datetime import datetime, timedelta, timezone

class PrayerCalculator:
    """
    Kalkulator Jadwal Sholat berbasis VSOP87D + rotasi bumi
    dengan konversi waktu UTC → WIB yang benar.
    """

    def __init__(self, vsop_path, lat, lon, timezone_offset):
        self.vsop_path = vsop_path
        self.lat = lat
        self.lon = lon
        self.tz = timezone_offset   # WIB = 7
        self.delta_t = 69.2 / 86400.0

    # ---------------- TIME HANDLING ----------------

    def local_to_jd(self, dt_local):
        dt_utc = dt_local - timedelta(hours=self.tz)
        return (dt_utc.timestamp() / 86400.0) + 2440587.5

    def jd_to_local_string(self, jd):
        dt_utc = datetime.fromtimestamp(
            (jd - 2440587.5) * 86400.0,
            tz=timezone.utc
        )
        dt_local = dt_utc + timedelta(hours=self.tz)
        return dt_local.strftime("%H:%M")

    # ---------------- SUN POSITION ----------------

    def get_sun_ra_dec(self, jd_ut):
        jd_tdb = jd_ut + self.delta_t
        sun = VSOP87D_Sun.get_sun_position(self.vsop_path, jd_tdb)
        return EarthRotation.ecliptic_to_equatorial(
            sun['lon'], sun['lat'], jd_tdb
        )

    # ---------------- CORE FALAK ----------------

    def solar_transit(self, jd_ut):
        ra, _ = self.get_sun_ra_dec(jd_ut)
        lst = (EarthRotation.get_gmst(jd_ut) + self.lon) % 360
        H = lst - math.degrees(ra)
        return jd_ut - H / 360.0

    def hour_angle(self, dec, alt):
        phi = math.radians(self.lat)
        alt = math.radians(alt)

        cosH = (
            math.sin(alt) - math.sin(phi) * math.sin(dec)
        ) / (
            math.cos(phi) * math.cos(dec)
        )

        if abs(cosH) > 1:
            return None

        return math.degrees(math.acos(cosH)) / 15.0

    # ---------------- MAIN CALC ----------------

    def calculate_daily_prayers(self, year, month, day):
        # Start langsung dari 12:00 WIB
        base_local = datetime(year, month, day, 12, 0, 0)
        jd_guess = self.local_to_jd(base_local)

        # Zuhur sejati
        zuhur_jd = self.solar_transit(jd_guess)

        # Deklinasi matahari
        _, dec = self.get_sun_ra_dec(zuhur_jd)

        # Ashar (Syafi'i)
        z = abs(math.radians(self.lat) - dec)
        ashar_alt = math.degrees(math.atan(1 / (1 + math.tan(z))))

        t_subuh   = self.hour_angle(dec, -20.0)
        t_maghrib = self.hour_angle(dec, -0.833)
        t_isya    = self.hour_angle(dec, -18.0)
        t_ashar   = self.hour_angle(dec, ashar_alt)

        return {
            "Subuh": self.jd_to_local_string(zuhur_jd - t_subuh/24),
            "Zuhur": self.jd_to_local_string(zuhur_jd),
            "Ashar": self.jd_to_local_string(zuhur_jd + t_ashar/24),
            "Maghrib": self.jd_to_local_string(zuhur_jd + t_maghrib/24),
            "Isya": self.jd_to_local_string(zuhur_jd + t_isya/24),
        }



if __name__ == "__main__":
    # Demo and Testing
    data_file = "VSOP87D.ear"
    jd_now = julian_now()
    
    print(f"Calculating for Julian Date (UT1): {jd_now:.6f}")
    
    # Sun Position
    sun = VSOP87D_Sun.get_sun_position(data_file, jd_now + 69.2/86400.0)
    print("\n--- GEOCENTRIC SUN POSITION ---")
    print(f"Longitude (Date): {sun['lon_deg']:.6f}°")
    
    # Prayer Times Example
    calc = PrayerCalculator(data_file, -6.2, 106.8, 7)
    times = calc.calculate_daily_prayers(2026, 2, 12)
    print("\n--- PRAYER TIMES (JAKARTA) ---")
    print(times)
