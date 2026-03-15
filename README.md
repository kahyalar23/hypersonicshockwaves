# 🚀 Hipersonik Şok Dalgaları — Hypersonic Shock Waves

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Mach 5+ hızlarda araç etrafındaki şok yapısı** — NASA Artemis programıyla doğrudan ilgili bir makine mühendisliği projesi.

---

## İçindekiler / Table of Contents

1. [Proje Hakkında / About](#proje-hakkında)
2. [Fizik / Physics](#fizik)
3. [Kurulum / Installation](#kurulum)
4. [Kullanım / Usage](#kullanım)
5. [Örnekler / Examples](#örnekler)
6. [Testler / Tests](#testler)
7. [Çıktılar / Output Figures](#çıktılar)
8. [Proje Yapısı / Project Structure](#proje-yapısı)
9. [Referanslar / References](#referanslar)

---

## Proje Hakkında

Bu proje, Mach 5 ve üzeri hızlarda (hipersonik rejim) bir araç etrafında oluşan
**şok dalgalarını** analiz eder ve görselleştirir.  NASA'nın **Artemis** programında
kullanılan **Orion kapsülü**, atmosfere yeniden giriş sırasında Mach ~32 hızına ulaşmakta
ve bu esnada devasa termodinamik yüklerle karşılaşmaktadır.

**Bu proje ile yapabilecekleriniz:**

- Normal şok dalgası (Rankine-Hugoniot) hesapları
- Eğik şok dalgası (θ–β–M ilişkisi)
- Prandtl-Meyer genişleme dalgaları
- Künt cisim önündeki detached bow shock şekli (Billig korelasyonu)
- Durma sıcaklığı ve basınç hesapları
- NASA Artemis yeniden giriş yörüngesi boyunca aerotermik analiz
- Mach 5–25 arası parametrik süpürme çalışmaları

---

## Fizik

### Normal Şok — Rankine-Hugoniot İlişkileri

| Büyüklük | Formül |
|---|---|
| M₂ | √[ ((γ−1)M₁² + 2) / (2γM₁² − (γ−1)) ] |
| p₂/p₁ | (2γM₁² − (γ−1)) / (γ+1) |
| ρ₂/ρ₁ | (γ+1)M₁² / ((γ−1)M₁² + 2) |
| T₂/T₁ | (p₂/p₁) / (ρ₂/ρ₁) |

### Eğik Şok — θ–β–M Bağıntısı

```
tan θ = 2 cot β · [M₁² sin²β − 1] / [M₁²(γ + cos 2β) + 2]
```

### Prandtl-Meyer Genişleme

```
ν(M) = √((γ+1)/(γ−1)) · arctan(√((γ−1)/(γ+1)·(M²−1))) − arctan(√(M²−1))
ν₂ = ν₁ + Δθ
```

### Bow Shock Standoff Mesafesi (Billig 1967)

```
δ/R = 0.386 · exp(4.67 / M∞²)
```

Detaylı teori için bkz. [`docs/theory.md`](docs/theory.md)

---

## Kurulum

```bash
# Repoyu klonlayın
git clone https://github.com/kahyalar23/hypersonicshockwavestest.git
cd hypersonicshockwavestest

# Bağımlılıkları yükleyin
pip install -r requirements.txt
```

**Gereksinimler:**

| Paket | Versiyon |
|---|---|
| numpy | ≥ 1.24 |
| scipy | ≥ 1.11 |
| matplotlib | ≥ 3.7 |

---

## Kullanım

### Python API

```python
from src.shock_relations import NormalShock, ObliqueShock
from src.prandtl_meyer import PrandtlMeyer
from src.bow_shock import BowShock

# Normal şok — Artemis doruk ısınma noktası (Mach 7.6)
ns = NormalShock(M1=7.6)
print(ns)
# NormalShock(M1=7.6, gamma=1.4)
#   M2        = 0.3963
#   p2/p1     = 66.6429
#   rho2/rho1 = 5.5628
#   T2/T1     = 11.9792
#   p02/p01   = 0.000210

# Eğik şok — kama yarı açısı θ = 20°
os_ = ObliqueShock(M1=10, theta=20)
print(f"Şok açısı β = {os_.beta:.2f}°")
print(f"p₂/p₁      = {os_.p_ratio:.2f}")

# Prandtl-Meyer genişlemesi
pm = PrandtlMeyer(M1=5, delta_theta=10)
print(f"Çıkış Mach  = {pm.M2:.3f}")

# Bow shock — Orion kapsülü (R_nose = 6 m)
bs = BowShock(M1=7.6, nose_radius=6.0)
print(bs)
```

---

## Örnekler

### 1. NASA Artemis Yeniden Giriş Analizi

```bash
python examples/artemis_reentry.py
```

Üç kritik yörünge noktasında (Mach 7.6, 10, 25) normal şok, bow shock ve
durma koşullarını hesaplar; çıktı grafikleri `output/` klasörüne kaydeder.

### 2. Parametrik Çalışma (Mach 5–25)

```bash
python examples/parametric_study.py
```

Normal şok özellikleri, eğik şok basınç oranı ısı haritası ve
bow shock standoff mesafesini Mach sayısının fonksiyonu olarak çizer.

### 3. Mach Konisi Görselleştirmesi

```bash
python examples/mach_cone.py
```

Mach 5, 10, 15, 20 için Mach konilerini ve bow shock şekillerini görselleştirir.

---

## Testler

```bash
pip install pytest
pytest tests/ -v
```

Toplam 40+ birim testi; Anderson'ın tablolarından alınan referans değerlerle
doğrulanmıştır.

---

## Çıktılar

Tüm örnekler çalıştırıldığında `output/` klasörüne aşağıdaki grafikler kaydedilir:

| Dosya | İçerik |
|---|---|
| `normal_shock_properties.png` | p, ρ, T, p₀ oranları vs. Mach |
| `artemis_bow_shocks.png` | Orion kapsülü — 3 yörünge noktasında bow shock |
| `artemis_summary.png` | Durma sıcaklığı ve basıncı özet grafiği |
| `stagnation_heating.png` | T₀ ve Taw vs. Mach (malzeme erime noktaları ile) |
| `parametric_normal_shock.png` | Normal şok parametrik analizi |
| `parametric_oblique_shock.png` | Eğik şok — 3 farklı θ için |
| `bow_shock_standoff.png` | Billig korelasyonu — δ/R vs. Mach |
| `theta_beta_mach.png` | Klasik θ–β–M diyagramı |
| `oblique_shock_heatmap.png` | p₂/p₁ ısı haritası (M, θ uzayı) |
| `mach_cone.png` | Mach konileri |
| `bow_shock_M7_6.png` | Artemis Mach 7.6 bow shock detayı |
| `bow_shock_M25.png` | Derin hipersonik Mach 25 bow shock |

---

## Proje Yapısı

```
hypersonicshockwavestest/
├── README.md
├── requirements.txt
├── docs/
│   └── theory.md              # Fizik referansı ve denklemler
├── src/
│   ├── __init__.py
│   ├── shock_relations.py     # Normal & eğik şok (Rankine-Hugoniot)
│   ├── prandtl_meyer.py       # Prandtl-Meyer genişleme dalgaları
│   ├── bow_shock.py           # Bow shock geometrisi (Billig korelasyonu)
│   └── visualization.py       # Görselleştirme araçları
├── examples/
│   ├── artemis_reentry.py     # NASA Artemis yeniden giriş vaka çalışması
│   ├── parametric_study.py    # Mach 5–25 parametrik analiz
│   └── mach_cone.py           # Mach konisi görselleştirmesi
├── tests/
│   ├── test_shock_relations.py
│   ├── test_prandtl_meyer.py
│   └── test_bow_shock.py
└── output/                    # Oluşturulan grafikler (git-ignored)
```

---

## NASA Artemis Bağlantısı

| Faz | Mach | Yükseklik | Yaklaşık T₀ |
|---|---|---|---|
| Giriş arayüzü | ~32 | 120 km | > 10 000 K |
| Doruk ısınma | ~10 | 75 km | ~3 200 K |
| Doruk dinamik basınç | ~7.6 | 65 km | ~1 800 K |
| Parashüt açılımı | < 0.8 | ~8 km | — |

Bu sıcaklıklar, çelik (1811 K) ve çoğu mühendislik alaşımını eriterek Orion'un
ablasyon korumalı ısı kalkanının neden gerekli olduğunu açıkça ortaya koymaktadır.

---

## Referanslar

1. Anderson, J. D. (2003). *Modern Compressible Flow*, 3rd ed. McGraw-Hill.
2. Anderson, J. D. (2006). *Hypersonic and High Temperature Gas Dynamics*, 2nd ed. AIAA.
3. Billig, F. S. (1967). Shock-wave shapes around spherical- and cylindrical-nosed bodies. *Journal of Spacecraft and Rockets*, 4(6), 822–823.
4. NACA TN 1428 — Equations, tables and charts for compressible flow (1953).
5. NASA SP-8077 — Aerodynamic Design Data Book.

---

## Lisans

MIT © 2024