# Hilal Tracker Falakiyah - Hisab Awal Bulan

Program interaktif berbasis CLI untuk melakukan hisab awal bulan Hijriah dengan akurasi tinggi menggunakan algoritma **VSOP87D** (Matahari) dan **ELP2000-82** (Bulan). Dirancang untuk memberikan data astronomis yang presisi bagi pengamat falak di Indonesia.

## 🚀 Fitur Utama

- **Akurasi Sub-Detik**: 
    - Menggunakan polinomial presisi tinggi dan koreksi **Delta-T** (Espenak & Meeus).
    - **Apparent Motion**: Memperhitungkan koreksi **Nutasi** dan **Aberasi** untuk akurasi waktu Ijtima dan posisi benda langit yang lebih tepat.
- **Hisab Seluruh Kota/Kab Indonesia**: 
    - Rekapitulasi otomatis untuk seluruh Kota dan Kabupaten di Indonesia sekaligus (Opsi 0).
    - Tampilan tabel yang rapi dengan fitur *text wrapping* untuk nama lokasi yang panjang.
    - Menampilkan kesimpulan jatuh tanggal 1 lengkap dengan **Nama Hari**.
- **Visualisasi Premium ASCII**:
    - **Twilight Gradient**: Menggunakan ANSI 256-color untuk simulasi gradasi langit senja.
    - **Side-by-Side View**: Menampilkan perbandingan posisi hilal antara hari ke-29 dan H+1 secara berdampingan.
    - **Dynamic Scaling**: Skala vertikal otomatis (auto-zoom) menyesuaikan tinggi hilal.
- **Kriteria MABIMS**: Deteksi otomatis berdasarkan kriteria terbaru (Alt >= 3°, Elong >= 6.4°).
- **Zona Waktu Administratif**: Mendukung pembagian WIB, WITA, dan WIT secara otomatis berdasarkan lokasi.

## 🛠️ Cara Penggunaan

1. Jalankan script utama:
   ```bash
   python hisab_awal_bulan.py
   ```
2. **Pilih Mode**:
   - Pilih nomor provinsi dan kota untuk detail spesifik lokasi.
   - Pilih `0` untuk **Hisab Seluruh Kota/Kab Indonesia**.
3. Masukkan Tahun Hijriyah yang ingin dihitung (e.g., 1447).
4. Navigasi:
   - Tekan nomor bulan (1-12) untuk melihat **Visualisasi ASCII Hilal**.
   - Tekan `Enter` untuk kembali ke menu utama.

## 📺 Demo & Tutorial

[![Tutorial Video](thumb.png)](https://github.com/user-attachments/assets/4cf8dc6e-562b-4e36-bdc5-97ce7af620e5)

## 📂 Struktur Project

- `hisab_awal_bulan.py`: Script utama (Logic UI, Visualisasi, dan Flow).
- `matahari.py`: Engine posisi Matahari (VSOP87D).
- `bulan.py`: Engine posisi Bulan (ELP2000-82).
- `cities_indonesia_processed.json`: Database koordinat dan elevasi lokasi di Indonesia.

---
**Kreator**: Tri Setyo Pamungkas  
**Engine**: VSOP87D & ELP2000-82  
**Website**: [hisab.tawheed.my.id](https://hisab.tawheed.my.id)
