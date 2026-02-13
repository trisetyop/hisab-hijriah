# Hilal Tracker Falakiyah - Hisab Awal Bulan

Program interaktif untuk melakukan hisab awal bulan Hijriah dengan akurasi tinggi menggunakan algoritma **VSOP87** untuk Matahari dan **ELP2000-82** untuk Bulan.

## Fitur Utama

- **Akurasi Tinggi**: Menggunakan polinomial presisi tinggi dan koreksi Delta-T (Espenak & Meeus). Seluruh kalkulasi diproses secara penuh **tanpa trunkasi** data.
- **Interaktif**: Pilih lokasi berdasarkan Provinsi dan Kota di seluruh Indonesia.
- **Zona Waktu Politik**: Mendukung pembagian WIB, WITA, dan WIT secara administratif.
- **Visualisasi ASCII**: Menampilkan posisi Hilal terhadap ufuk secara visual di terminal.
- **Kriteria MABIMS**: Deteksi otomatis awal bulan berdasarkan kriteria baru MABIMS (Alt >= 3°, Elong >= 6.4°).
- **Dukungan Sejarah**: Mendukung kalender Julian (sebelum 1582) dan Gregorian.

## Cara Penggunaan

1. Jalankan script utama:
   ```bash
   python hisab_awal_bulan.py
   ```
2. Pilih Provinsi dan Kota.
3. Masukkan tahun Hijriah yang ingin dihitung.
4. Lihat tabel hasil hisab selama 12 bulan.
5. Pilih nomor bulan (1-12) untuk melihat visualisasi detail Hilal.

## Tutorial Video

Tonton panduan penggunaan lengkap di sini:
[![Tutorial Video](thumb.png)](https://github.com/user-attachments/assets/4cf8dc6e-562b-4e36-bdc5-97ce7af620e5)


## Struktur File Penting

- `hisab_awal_bulan.py`: Script utama program.
- `matahari.py`: Engine posisi Matahari (VSOP87).
- `bulan.py`: Engine posisi Bulan (ELP2000-82).
- `provinces.json`: Data provinsi Indonesia.
- `cities_indonesia_processed.json`: Data koordinat kota-kota di Indonesia.

---
**Kreator**: Tri Setyo Pamungkas  
**Website**: [hisab.tawheed.my.id](https://hisab.tawheed.my.id)
