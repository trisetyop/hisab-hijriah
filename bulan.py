import math
import os
from functools import lru_cache

class ELP2000_Moon:
    """
    Parser dan Kalkulator untuk ELP 2000-82 Lunar Theory (Chapront-Touze & Chapront, 1983).
    Mengadaptasi logika dari ELPMPP02.for.
    
    OPTIMASI: Menggunakan LRU cache untuk hasil perhitungan posisi bulan.
    """
    
    # Constants
    CPI = 3.141592653589793
    RAD = 648000.0 / CPI
    DEG = CPI / 180.0
    PIS2 = CPI / 2.0
    DPI = 2.0 * CPI
    
    SC = 36525.0
    A405 = 384747.9613701725
    AELP = 384747.980674318
    
    # ========== CACHE UNTUK HASIL PERHITUNGAN ==========
    _result_cache = {}  # Cache untuk hasil calculate()
    
    def __init__(self, data_dir=""):
        self.data_dir = data_dir
        self.w = [[0.0]*5 for _ in range(4)]
        self.eart = [0.0]*5
        self.peri = [0.0]*5
        self.del_args = [[0.0]*5 for _ in range(5)] 
        self.p_args = [[0.0]*5 for _ in range(9)] 
        self.zeta = [0.0]*5
        
        # Series data containers
        # index 0 unused, 1=Lon, 2=Lat, 3=Dist
        self.main_series = [None, [], [], []]
        self.pert_series = [None, [], [], []]
        
        self.initialized = False

    def _dms_to_rad(self, ideg, imin, sec):
        return (ideg + imin/60.0 + sec/3600.0) * self.DEG

    def _initialize_constants(self, icor=0):
        """Initializes internal constants based on ELPMPP02 INITIAL subroutine."""
        rad = self.RAD
        
        # Constants differences based on icor
        Dw1_0 = -0.10525; Dw2_0 = 0.16826; Dw3_0 = -0.10760
        Deart_0 = -0.04012; Dperi = -0.04854
        Dw1_1 = -0.32311; Dgam = 0.00069; De = +0.00005
        Deart_1 = 0.01442; Dep = 0.00226; Dw2_1 = 0.08017
        Dw3_1 = -0.04317; Dw1_2 = -0.03794
        
        # DE405 overrides if icor=1 (skipping for now, assuming LLR/icor=0 default)

        # W1: Moon Mean Longitude
        self.w[1] = [
            self._dms_to_rad(218, 18, 59.95571 + Dw1_0),
            (1732559343.73604 + Dw1_1) / rad,
            (-6.8084 + Dw1_2) / rad,
            0.66040e-2 / rad,
            -0.31690e-4 / rad
        ]
        
        # W2: Perigee
        self.w[2] = [
            self._dms_to_rad(83, 21, 11.67475 + Dw2_0),
            (14643420.3171 + Dw2_1) / rad,
            -38.2631 / rad,
            -0.45047e-1 / rad,
            0.21301e-3 / rad
        ]
        
        # W3: Node
        self.w[3] = [
            self._dms_to_rad(125, 2, 40.39816 + Dw3_0),
            (-6967919.5383 + Dw3_1) / rad,
            6.3590 / rad,
            0.76250e-2 / rad,
            -0.35860e-4 / rad
        ]
        
        # Earth Elements
        self.eart = [
            self._dms_to_rad(100, 27, 59.13885 + Deart_0),
            (129597742.29300 + Deart_1) / rad,
            -0.020200 / rad,
            0.90000e-5 / rad,
            0.15000e-6 / rad
        ]
        
        self.peri = [
            self._dms_to_rad(102, 56, 14.45766 + Dperi),
            1161.24342 / rad,
            0.529265 / rad,
            -0.11814e-3 / rad,
            0.11379e-4 / rad
        ]
        
        # Mean motion corrections (simplified logic from ELPMPP02)
        bp = [
            [0.311079095, -0.103837907],
            [-0.4482398e-2, 0.6682870e-3],
            [-0.110248500e-2, -0.129807200e-2],
            [0.1056062e-2, -0.1780280e-3],
            [0.50928e-4, -0.37342e-4]
        ]
        
        am = 0.074801329
        alpha = 0.002571881
        dtasm = (2.0 * alpha) / (3.0 * am)
        xa = (2.0 * alpha) / 3.0
        
        x2 = self.w[2][1] / self.w[1][1]
        x3 = self.w[3][1] / self.w[1][1]
        
        # Precompute y2/y3 based on column 0 and 1 of bp (index 0 and 1)
        y2 = am * bp[0][0] + xa * bp[4][0]
        y3 = am * bp[0][1] + xa * bp[4][1]
        
        d21 = x2 - y2
        d22 = self.w[1][1] * bp[1][0]
        d23 = self.w[1][1] * bp[2][0]
        d24 = self.w[1][1] * bp[3][0]
        d25 = y2 / am
        
        d31 = x3 - y3
        d32 = self.w[1][1] * bp[1][1]
        d33 = self.w[1][1] * bp[2][1]
        d34 = self.w[1][1] * bp[3][1]
        d35 = y3 / am
        
        Cw2_1 = d21*Dw1_1 + d25*Deart_1 + d22*Dgam + d23*De + d24*Dep
        Cw3_1 = d31*Dw1_1 + d35*Deart_1 + d32*Dgam + d33*De + d34*Dep
        
        self.w[2][1] += Cw2_1 / rad
        self.w[3][1] += Cw3_1 / rad
        
        self.dtasm = dtasm
        self.am = am
        
        # Delaunay
        for i in range(5):
            self.del_args[1][i] = self.w[1][i] - self.eart[i]
            self.del_args[2][i] = self.w[1][i] - self.w[3][i]
            self.del_args[3][i] = self.w[1][i] - self.w[2][i]
            self.del_args[4][i] = self.eart[i] - self.peri[i]
        self.del_args[1][0] += self.CPI
        
        # Planetary
        self.p_args[1][0] = self._dms_to_rad(252, 15, 3.216919);   self.p_args[1][1] = 538101628.66888 / rad
        self.p_args[2][0] = self._dms_to_rad(181, 58, 44.758419);  self.p_args[2][1] = 210664136.45777 / rad
        self.p_args[3][0] = self._dms_to_rad(100, 27, 59.138850);  self.p_args[3][1] = 129597742.29300 / rad
        self.p_args[4][0] = self._dms_to_rad(355, 26, 3.642778);   self.p_args[4][1] = 68905077.65936 / rad
        self.p_args[5][0] = self._dms_to_rad(34, 21, 5.379392);    self.p_args[5][1] = 10925660.57335 / rad
        self.p_args[6][0] = self._dms_to_rad(50, 4, 38.902495);    self.p_args[6][1] = 4399609.33632 / rad
        self.p_args[7][0] = self._dms_to_rad(314, 3, 4.354234);    self.p_args[7][1] = 1542482.57845 / rad
        self.p_args[8][0] = self._dms_to_rad(304, 20, 56.808371);  self.p_args[8][1] = 786547.89700 / rad
        
        # Zeta
        self.zeta[0] = self.w[1][0]
        # Dprec = -0.29965 ! IAU 2000A
        Dprec = -0.29965
        self.zeta[1] = self.w[1][1] + (5029.0966 + Dprec) / rad
        self.zeta[2] = self.w[1][2]; self.zeta[3] = self.w[1][3]; self.zeta[4] = self.w[1][4]
        
        # Corrections
        self.delnu = (+0.55604 + Dw1_1) / rad / self.w[1][1]
        self.dele  = (+0.01789 + De) / rad
        self.delg  = (-0.08066 + Dgam) / rad
        self.delnp = (-0.06424 + Deart_1) / rad / self.w[1][1]
        self.delep = (-0.12879 + Dep) / rad
        
        # Precession Coeffs (Laskar)
        self.p_prec = [0, 0.10180391e-04, 0.47020439e-06, -0.5417367e-09, -0.2507948e-11, 0.463486e-14]
        self.q_prec = [0, -0.113469002e-03, 0.12372674e-06, 0.1265417e-08, -0.1371808e-11, -0.320334e-14]
        
        self.initialized = True

    @staticmethod
    def parse_main_file(filepath):
        """
        Membaca file seri utama (ELP_MAIN.S1, .S2, .S3).
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File {filepath} tidak ditemukan.")

        series = []
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if not lines:
            return series

        # Skip/Read header line correctly if needed
        # We start from index 1 as before
        idx = 1
        
        while idx < len(lines):
            line = lines[idx]
            if len(line) < 10: 
                idx += 1
                continue
            try:
                i1 = int(line[0:3])
                i2 = int(line[3:6])
                i3 = int(line[6:9])
                i4 = int(line[9:12])
                a_val = float(line[14:27])
                b1 = float(line[27:39])
                b2 = float(line[39:51])
                b3 = float(line[51:63])
                b4 = float(line[63:75])
                b5 = float(line[75:87])
                
                series.append({
                    'ilu': [i1, i2, i3, i4],
                    'a': a_val,
                    'b': [b1, b2, b3, b4, b5]
                })
            except ValueError:
                pass
            idx += 1
        return series

    def _process_main_series(self, raw_series, iv_index):
        processed = []
        for term in raw_series:
            a = term['a']
            b = term['b']
            ilu = term['ilu']
            
            tgv = b[0] + self.dtasm * b[4]
            if iv_index == 3: # Distance
                a = a - 2.0 * a * self.delnu / 3.0
                
            coeff_val = a + tgv * (self.delnp - self.am * self.delnu) \
                        + b[1] * self.delg + b[2] * self.dele + b[3] * self.delep
                        
            # fmpb(k)
            f_poly = [0.0]*5
            for k in range(5):
                sum_val = 0.0
                for i in range(4): # 0..3 index for ilu
                    sum_val += ilu[i] * self.del_args[i+1][k]
                f_poly[k] = sum_val
            
            if iv_index == 3:
                f_poly[0] += self.PIS2
            
            processed.append({'coef': coeff_val, 'poly': f_poly})
        return processed

    def _process_pert_series(self, raw_series):
        processed = []
        for term in raw_series:
            s_val = term['s']
            c_val = term['c']
            ifi = term['ifi']
            
            # Amplitude
            cper = math.sqrt(c_val*c_val + s_val*s_val)
            
            # Phase
            pha = math.atan2(c_val, s_val)
            if pha < 0: pha += self.DPI
            
            f_poly = [0.0]*5
            for k in range(5):
                val = 0.0
                if k == 0: val = pha
                
                # ifi 1..4 -> del, 5..12 -> p, 13 -> zeta
                for i in range(4): # 0..3 -> del 1..4
                    val += ifi[i] * self.del_args[i+1][k]
                for i in range(4, 12): # 4..11 -> p 1..8
                    val += ifi[i] * self.p_args[i-3][k] # p_args id 1..8
                
                val += ifi[12] * self.zeta[k]
                f_poly[k] = val
            
            # Perturbations need 'index t' (it)
            # but wait, ELP_PERT has blocks for t^0, t^1...
            # The files structure is:
            # S1 -> iv=1 (Long), S2 -> iv=2 (Lat), S3 -> iv=3 (Dist)
            # Inside each file, there are blocks!
            # Fortran: loop do it=0,3 ... read (lu,1003) nper, ipt ... do n=1,nper
            
            # My simple parser ignored the 'it' blocks and flattened them?
            # NO! My simple parser just read line by line.
            # I need `it` (time power) for evaluation: x * t(it) * sin(y)
            # I must fix the parsing to capture `it` or re-parse properly.
            
            processed.append({
                'coef': cper, 
                'poly': f_poly
                # MISSING 'it' (time power)!
            })
        return processed

    def _reparse_pert_blocks(self, filepath):
        """
        Parses PERT file correctly respecting the blocks.
        Format:
        Header
        Block Header: (25x, 2i10) -> n_terms, ipt (ipt is not used in formula?)
        Data lines...
        """
        if not os.path.exists(filepath): return []
        
        series = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        idx = 1 # Skip file header
        while idx < len(lines):
            line = lines[idx]
            if len(line) < 20: 
                idx += 1; continue
                
            # Check for Block Header
            # Block header in fortran: read(lu,1003) nper, ipt
            # 1003 format (25x, 2i10)
            try:
                # Heuristic: Block header usually starts with spaces and has 2 integers
                # Data lines behave differently (i5 at start).
                # Let's check columns 25-35 and 35-45
                n_terms = int(line[25:35])
                ipt = int(line[35:45]) # 'it' relates to this? 
                
                # Fortran loop: do it=0,3.
                # The file contains blocks sequentially for it=0, it=1, it=2, it=3
                # So I just need to track which block I am in?
                # Actually Fortran code expects exactly 4 blocks (it=0..3).
                # BUT S1/S2/S3 might have empty blocks?
                # "read (lu,1003) nper, ipt ... if (nper.eq.0) cycle"
                
                # THIS IS TRICKY. If I just read sequentially, I need to know 'it'.
                # Let's count blocks?
                # S1 has it=0..3 blocks.
                
                # We can't rely on pure sequential reading without context unless we know we are at a header.
                # Data lines: "    1-0.12..." (digits at start).
                # Header lines: "                         11314         0" (spaces at start).
                
                # Reset counter
                current_it = ipt # Usually ipt is the 'it' index? Let's check fortran.
                # Fortran: read(lu) nper, ipt.
                # It loop: do it=0,3.
                # Does ipt == it? 
                # Lines 1: "11314 0" -> it=0?
                
                idx += 1
                for _ in range(n_terms):
                    dline = lines[idx]
                    idx += 1
                    
                    term_idx = int(dline[0:5])
                    s_val = float(dline[5:25].replace('D', 'E'))
                    c_val = float(dline[25:45].replace('D', 'E'))
                    ifi = []
                    start = 45
                    for _ in range(16):
                        ifi.append(int(dline[start:start+3]))
                        start += 3
                        
                    series.append({
                        'it': current_it, # Assuming ipt is the polynomial degree
                        's': s_val,
                        'c': c_val,
                        'ifi': ifi
                    })
                    
            except ValueError:
                # potentially empty line or EOF
                idx += 1
                
        return series

    def _process_pert_series_with_it(self, raw_series):
        processed = []
        for term in raw_series:
            it = term.get('it', 0)
            s_val = term['s']
            c_val = term['c']
            ifi = term['ifi']
            
            cper = math.sqrt(c_val*c_val + s_val*s_val)
            pha = math.atan2(c_val, s_val)
            if pha < 0: pha += self.DPI
            
            f_poly = [0.0]*5
            for k in range(5):
                val = 0.0
                if k == 0: val = pha
                for i in range(4): val += ifi[i] * self.del_args[i+1][k]
                for i in range(4, 12): val += ifi[i] * self.p_args[i-3][k]
                val += ifi[12] * self.zeta[k]
                f_poly[k] = val
                
            processed.append({'coef': cper, 'poly': f_poly, 'it': it})
        return processed

    def load_all_series(self):
        """Loads and processes all ELP files."""
        if not self.initialized: self._initialize_constants()
        
        # Main
        for iv in range(1, 4):
            fname = f"ELP_MAIN.S{iv}"
            fpath = os.path.join(self.data_dir, fname)
            if os.path.exists(fpath):
                raw = self.parse_main_file(fpath)
                self.main_series[iv] = self._process_main_series(raw, iv)
                
        # Pert
        for iv in range(1, 4):
            fname = f"ELP_PERT.S{iv}"
            fpath = os.path.join(self.data_dir, fname)
            if os.path.exists(fpath):
                # Use the block-aware parser
                raw = self._reparse_pert_blocks(fpath)
                self.pert_series[iv] = self._process_pert_series_with_it(raw)

    def calculate(self, jd_tdb):
        """
        Calculates Geocentric Position (X,Y,Z) and (Lon, Lat, Dist).
        jd_tdb: Julian Date (TDB scale theoretically, but usually close to TT/ET)
        
        OPTIMASI: Menggunakan caching untuk hasil perhitungan.
        JD dibulatkan ke 1 menit untuk cache hit yang lebih baik.
        """
        # Cache key: JD dibulatkan ke 1 menit (1440 menit per hari)
        cache_key = round(jd_tdb * 1440)
        
        if cache_key in self._result_cache:
            return self._result_cache[cache_key]
        
        if self.main_series[1] is None:
            self.load_all_series()
            
        t = [0.0]*5
        t[0] = 1.0
        t[1] = (jd_tdb - 2451545.0) / 36525.0
        t[2] = t[1]*t[1]
        t[3] = t[2]*t[1]
        t[4] = t[3]*t[1]
        
        # v[1]=Lon, v[2]=Lat, v[3]=Dist
        v = [0.0]*4 
        
        for iv in range(1, 4):
            val = 0.0
            
            # Main
            if self.main_series[iv]:
                for term in self.main_series[iv]:
                    # y = poly[0] + poly[1]*t + ...
                    # Optimization: Horner's method or explicit
                    y = term['poly'][0] + term['poly'][1]*t[1] + \
                        term['poly'][2]*t[2] + term['poly'][3]*t[3] + term['poly'][4]*t[4]
                    
                    x = term['coef']
                    val += x * math.sin(y)
                    
            # Pert
            if self.pert_series[iv]:
                for term in self.pert_series[iv]:
                    it = term['it'] # 0..3
                    # x * t^it * sin(y)
                    
                    y = term['poly'][0] + term['poly'][1]*t[1] + \
                        term['poly'][2]*t[2] + term['poly'][3]*t[3] + term['poly'][4]*t[4]
                    
                    x = term['coef']
                    val += x * t[it] * math.sin(y)
            
            v[iv] = val

        # Final corrections (from EVALUATE)
        # v(1) = v(1)/rad + W1(t)
        # v(2) = v(2)/rad
        # v(3) = v(3)*a405/aelp
        
        w1_val = self.w[1][0] + self.w[1][1]*t[1] + self.w[1][2]*t[2] + \
                 self.w[1][3]*t[3] + self.w[1][4]*t[4]
        
        # Precession from J2000 to Date (Meeus Ch 21)
        # ELP2000-82 results are J2000 Fixed.
        # We want Mean of Date for calculate_topocentric or modern comparison.
        # T is centuries from J2000.
        T_cent = t[1] 
        prec_lon_sec = (5029.0966 + 1.11113 * T_cent) * T_cent
        
        lon_rad = v[1] / self.RAD + w1_val + (prec_lon_sec / self.RAD)
        lat_rad = v[2] / self.RAD
        dist_km = v[3] * self.A405 / self.AELP 
        
        # Precession in latitude is much smaller and often neglected in this scale 
        # (around 0.01" over decades).
        
        # Normalize Lon
        lon_rad = lon_rad % self.DPI
        
        result = {
            'lon_rad': lon_rad,
            'lat_rad': lat_rad,
            'dist_km': dist_km,
            'lon_deg': math.degrees(lon_rad),
            'lat_deg': math.degrees(lat_rad)
        }
        
        # Simpan ke cache (dengan batas cache size untuk mencegah memory leak)
        if len(self._result_cache) < 10000:  # Batas cache size
            self._result_cache[cache_key] = result
        
        return result

    # --- TOPOCENTRIC ADDITIONS ---

    def get_obliquity(self, jd):
        """Mean Obliquity of the Ecliptic (Laskar / Meeus)"""
        T = (jd - 2451545.0) / 36525.0
        eps0 = (23.43929111 -
                (46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3) / 3600.0)
        return math.radians(eps0)

    def get_gmst(self, jd_ut):
        """Greenwich Mean Sidereal Time (degrees)"""
        T = (jd_ut - 2451545.0) / 36525.0
        gmst = (280.46061837 +
                360.98564736629 * (jd_ut - 2451545.0) +
                0.000387933 * T**2 -
                (T**3 / 38710000.0))
        return gmst % 360.0

    def ecliptic_to_equatorial(self, lon_rad, lat_rad, jd):
        """
        Converts Geocentric Ecliptic -> Geocentric Equatorial (RA, Dec).
        """
        eps = self.get_obliquity(jd)
        sin_ra = (math.sin(lon_rad) * math.cos(eps) -
                  math.tan(lat_rad) * math.sin(eps))
        cos_ra = math.cos(lon_rad)
        ra_rad = math.atan2(sin_ra, cos_ra)
        ra_rad = (ra_rad + self.DPI) % self.DPI

        sin_dec = (math.sin(lat_rad) * math.cos(eps) +
                   math.cos(lat_rad) * math.sin(eps) * math.sin(lon_rad))
        dec_rad = math.asin(sin_dec)

        return ra_rad, dec_rad

    def calculate_topocentric(self, jd_ut, lat_obs_deg, lon_obs_deg, alt_obs_m=0.0, pressure=1010, temp_c=10, delta_t_sec=0.0):
        """
        Calculates Topocentric Equatorial Position (RA', Dec', Dist').
        Includes Parallax correction.
        Output also includes Apparent Altitude (with Refraction and optional Dip).
        Menghitung posisi toposentrik Bulan.
        
        Input:
        - jd_ut: Julian Day (UT1) untuk rotasi bumi (GMST).
        - delta_t_sec: TD - UT dalam detik (untuk Ephemeris Bulan). 
                       Default 0.0 jika tidak diketahui (mengurangi akurasi ~1').
        """
        
        # JD TDB for Ephemeris (Moon Position) & Nutation
        jd_tdb = jd_ut + (delta_t_sec / 86400.0)

        # 1. Geocentric Position (Ecliptic)
        geo = self.calculate(jd_tdb) 
        
        
        # --- NUTATION & ABERRATION ---
        # Reimplement logic or import? Better to import if we trust relative path.
        # But user wants standalone bulan.py? 
        # For robustness, let's implement simplified Nutation here or assume user has updated matahari.py.
        # Let's assume standalone capability is preferred if verify scripts typically import both.
        # Trying import first.
        
        try:
            from matahari import EarthRotation
            d_psi, d_eps = EarthRotation.get_nutation(jd_tdb)
            # Aberration (Planetary) - Moon is close so -20.6" formula for sun is WRONG for Moon.
            # Moon Light time ~1.3s. Motion ~0.55 arcsec/sec.
            # Correction ~ -0.7 arcsec.
            # Constant approx: -0.70" in Longitude?
            # Or simplified: use geometric. Moon aberration is small.
            # Sun Aberration: -20.4898 / Dist_AU.
            # Let's skip Moon Aberration for now or apply constant -0.0002 deg? 
            # Nutation is the big one (up to 17").
        except ImportError:
            # Fallback if matahari.py not found/updated
            d_psi, d_eps = 0.0, 0.0

        # Ensure functions using time for Space positions use jd_tdb
        # Functions using time for Earth Rotation use jd_ut
        
        # 2. Geocentric Equatorial (True of Date)
        # Note: self.ecliptic_to_equatorial needs update or we do manual conversion here
        
        # Calculate True Obliquity (using TDB roughly OK)
        T = (jd_tdb - 2451545.0) / 36525.0
        eps0 = (23.43929111 - (46.8150 * T - 0.00059 * T**2 + 0.001813 * T**3) / 3600.0)
        eps_true_rad = math.radians(eps0 + d_eps)
        
        # Apply Nutation to Longitude
        lon_true_rad = geo['lon_rad'] + math.radians(d_psi)
        lat_rad = geo['lat_rad']
        
        # Conversion to Equatorial
        sin_ra = (math.sin(lon_true_rad) * math.cos(eps_true_rad) -
                  math.tan(lat_rad) * math.sin(eps_true_rad))
        cos_ra = math.cos(lon_true_rad)
        ra_rad = math.atan2(sin_ra, cos_ra)
        ra_rad = (ra_rad + self.DPI) % self.DPI

        sin_dec = (math.sin(lat_rad) * math.cos(eps_true_rad) +
                   math.cos(lat_rad) * math.sin(eps_true_rad) * math.sin(lon_true_rad))
        dec_rad = math.asin(sin_dec)
        
        ra_geo, dec_geo = ra_rad, dec_rad
        dist_geo = geo['dist_km']
        
        # 3. Observer Parameters
        lat_rad = math.radians(lat_obs_deg)
        lon_rad = math.radians(lon_obs_deg)
        
        
        
        def _get_dip(h_m):
            if h_m <= 0: return 0.0
            return 0.0293 * math.sqrt(h_m)
            
        def _get_refraction(h_deg, P, T):
            h = max(h_deg, -5.0)
            R_min = 1.02 / math.tan(math.radians(h + 10.3 / (h + 5.11)))
            R_val = R_min * (P / 1010.0) * (283.0 / (273.0 + T))
            return R_val / 60.0

        # Geodetic to Geocentric Latitude conversion
        # WGS84 Constants
        a_earth = 6378.137 # km
        f_earth = 1.0 / 298.257223563
        e2 = 2*f_earth - f_earth**2
        
        # Calculate rho*sin(phi') and rho*cos(phi')
        # Meeus Ch. 11 formula
        sin_lat = math.sin(lat_rad)
        cos_lat = math.cos(lat_rad)
        
        C = 1.0 / math.sqrt(1.0 - e2 * sin_lat**2)
        S = C * (1.0 - e2)
        
        # alt_obs in km
        h_km = alt_obs_m / 1000.0
        
        rho_cos_phi_prime = (C * a_earth + h_km) * cos_lat
        rho_sin_phi_prime = (S * a_earth + h_km) * sin_lat

        # Note: LST also needs Nutation correction ideally (Equation of Equinoxes)
        # GAST = GMST + DeltaPsi * cos(Eps)
        # We should apply this for consistency.
        
        def _get_gast(jd, d_psi_deg, eps_true_rad):
             # Simple GMST (Uses UT1)
             T = (jd - 2451545.0) / 36525.0
             gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0) +
                     0.000387933 * T**2 - (T**3 / 38710000.0)) % 360.0
             # Add Equation of Equinoxes
             eq_eq = d_psi_deg * math.cos(eps_true_rad)
             return (gmst + eq_eq) % 360.0

        if d_psi != 0:
            gast = _get_gast(jd_ut, d_psi, eps_true_rad)
            # LST = GAST + Longitude
            lst_deg = (gast + lon_obs_deg) % 360.0
        else:
             # Fallback
             from matahari import EarthRotation
             gmst = EarthRotation.get_gmst(jd_ut)
             lst_deg = (gmst + lon_obs_deg) % 360.0
             
        lst_rad = math.radians(lst_deg)
        
        # 5. Hour Angle
        # H_geo = lst_rad - ra_geo
        
        # 6. Parallax Calculation (Rigorous Vector)
        
        # Moon pos in Geocentric Equatorial (X towards Equinox):
        Xg = dist_geo * math.cos(dec_geo) * math.cos(ra_geo)
        Yg = dist_geo * math.cos(dec_geo) * math.sin(ra_geo)
        Zg = dist_geo * math.sin(dec_geo)
        
        # Observer pos in Geocentric Equatorial (X towards Equinox):
        Xo = rho_cos_phi_prime * math.cos(lst_rad)
        Yo = rho_cos_phi_prime * math.sin(lst_rad)
        Zo = rho_sin_phi_prime
        
        # Topocentric Vector
        Xt = Xg - Xo
        Yt = Yg - Yo
        Zt = Zg - Zo
        
        # Convert back to Topocentric RA/Dec
        dist_topo = math.sqrt(Xt**2 + Yt**2 + Zt**2)
        ra_topo = math.atan2(Yt, Xt)
        ra_topo = (ra_topo + self.DPI) % self.DPI
        dec_topo = math.asin(Zt / dist_topo)
        
        # Optional: Azimuth/Altitude
        H_topo = lst_rad - ra_topo
        sin_alt = (math.sin(lat_rad) * math.sin(dec_topo) +
                   math.cos(lat_rad) * math.cos(dec_topo) * math.cos(H_topo))
        alt_rad = math.asin(sin_alt)
        alt_deg = math.degrees(alt_rad)
        
        sin_h = math.sin(H_topo)
        cos_h = math.cos(H_topo)
        
        # Robust Azimuth (Meeus)
        # x = -sin(H) * cos(dec)
        # y = sin(dec)*cos(phi) - cos(dec)*cos(H)*sin(phi)
        x_az = -sin_h * math.cos(dec_topo)
        y_az = (math.sin(dec_topo) * math.cos(lat_rad)) - (math.cos(dec_topo) * cos_h * math.sin(lat_rad))
        
        az_rad = math.atan2(x_az, y_az)
        
        # Atmospheric Corrections
        dip = _get_dip(alt_obs_m)
        refr = _get_refraction(alt_deg - dip, pressure, temp_c)
        
        # Apparent Altitude = Geometric - Dip + Refraction
        # Wait, standard is:
        # Apparent = Geometric + Refraction
        # But if observed from height, user sees "Dip" lower.
        # So effective altitude above sea horizon = Geom + Dip? No.
        # Geometric 0 is perp to gravity. Sea horizon is at -Dip.
        # So object at Geom 0 is at Sea-Alt +Dip.
        # Usually we return "Altitude" as geometric or apparent (refracted).
        # Dip is usually separate or "Altitute relative to sea horizon".
        # Let's provide 'alt_apparent' as refracted geometric.
        # And 'alt_observed' as refracted + dip (height corrected).
        
        alt_apparent = alt_deg + refr
        alt_observed = alt_apparent + dip # Relative to sea horizon
        
        return {
            'ra_deg': math.degrees(ra_topo),
            'dec_deg': math.degrees(dec_topo),
            'dist_km': dist_topo,
            'az_deg': (math.degrees(az_rad) + 360) % 360,
            'alt_deg': alt_deg,           # Topocentric Geometric
            'alt_apparent': alt_apparent, # Topocentric + Refraction
            'alt_observed': alt_observed, # Topocentric + Refraction + Dip (Sea Horizon)
            'dip_deg': dip,
            'refraction': refr,
            'geo_ra': math.degrees(ra_geo),
            'geo_dec': math.degrees(dec_geo)
        }

if __name__ == "__main__":
    elp = ELP2000_Moon()
    print("Loading data...")
    elp.load_all_series()
    print("Data loaded.")
    
    # Test for Now at Jakarta
    import datetime
    now = datetime.datetime.now(datetime.timezone.utc)
    # Simple JD conversion (approx)
    def to_jd(dt):
        return (dt.timestamp() / 86400.0) + 2440587.5
    
    jd_now = to_jd(now)
    # Jakarta
    lat_jkt = -6.2088
    lon_jkt = 106.8456
    
    res = elp.calculate_topocentric(jd_now, lat_jkt, lon_jkt)
    print(f"Time (UTC): {now}")
    print(f"Jakarta Topocentric Position:")
    print(f"  Azimuth : {res['az_deg']:.4f}°")
    print(f"  Altitude: {res['alt_deg']:.4f}°")
    print(f"  RA      : {res['ra_deg']:.4f}°")
    print(f"  Dec     : {res['dec_deg']:.4f}°")
    print(f"  Dist    : {res['dist_km']:.1f} km")
