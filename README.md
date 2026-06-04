# BLM3810 Biyoenformatiğe Giriş - Akciğer Adenokarsinomu (GSE68465) Projesi

Bu proje, Yıldız Teknik Üniversitesi Bilgisayar Mühendisliği Bölümü **BLM3810 Biyoenformatiğe Giriş** dersi kapsamında geliştirilmiştir. Proje, Akciğer Adenokarsinomu (Lung Adenocarcinoma) hastalarından alınan gen ifadesi verileri (NCBI GEO: GSE68465) üzerinde yapılan biyoenformatik analizleri (Aşama 1) ve hastaların histolojik tümör derecelerini (Well Differentiated vs. Poorly Differentiated) tahmin eden makine öğrenmesi modellerini (Aşama 2) içermektedir.

Sistem yükleme limitleri nedeniyle ham veri setleri (boyutları çok büyük olduğu için) temizlenmiştir. Ön işleme kodunu çalıştırdığınızda eksik olan ham veriler NCBI GEO üzerinden otomatik olarak tekrar indirilecektir.

---

## 📂 Proje Dizin Yapısı

*   📁 **`data/`**: Temizlenmiş ve ön işlemeden geçmiş RDS veri dosyalarını (`GSE68465.rds` ve `GSE68465_label.rds`) barındırır.
*   📁 **`plots/`**: Çalıştırılan scriptlerin ürettiği tüm görselleştirmeleri (Boxplot, Keman grafiği, Kaplan-Meier eğrileri, WGCNA grafikleri, GO zenginleştirme grafikleri, ROC eğrileri ve Hata matrisleri) barındırır.
*   📁 **`results/`**: Diferansiyel ifade tablolarını, modül gen listelerini, WGCNA hub genlerini, Cytoscape düğüm/kenar dosyalarını ve makine öğrenmesi performans karşılaştırma tablolarını barındırır.
*   📄 **`packages.R`**: Projede kullanılan tüm R ve Bioconductor kütüphanelerinin kurulumunu yapan R scripti.
*   📄 **`GSE68465_preprocessing.R`**: Veri setini GEO'dan indiren, temizleyen, prob-gen eşlemesi yapan ve klinik meta-veriyi düzenleyen ön işleme scripti.
*   📄 **`GSE68465_DGE.R`**: Histolojik dereceye göre (Well vs Poorly) diferansiyel ifade gösteren genleri `limma` ile analiz eden script.
*   📄 **`GSE68465_visualization.R`**: Volkan (Volcano) grafiği ile en anlamlı genlerin Wilcoxon p-değerli tekli/çoklu boxplot ve violin grafiklerini çizen script.
*   📄 **`GSE68465_survival.R`**: En anlamlı genler için Kaplan-Meier eğrilerini çıkaran ve klinik değişkenlerle çok değişkenli Cox regresyonu yapıp Forest plot çizen script.
*   📄 **`GSE68465_wgcna.R`**: Ağırlıklı Gen Eş-ifade Ağı Analizi (WGCNA) yaparak gen modüllerini çıkaran, hub genleri bulan ve Cytoscape ağ dosyalarını üreten script.
*   📄 **`GSE68465_GO.R`**: Histolojik dereceyle en çok ilişkili gen modülü için Gene Ontology (GO) zenginleştirme analizi yapan script.
*   📄 **`GSE68465_ML.R`**: Aşama 2 gereksinimlerini karşılayan, 3 farklı özellik seçimi düzeneği ve 3 makine öğrenmesi modeli (Random Forest, SVM, Lojistik Regresyon) ile 5-Fold çapraz doğrulama yapan script.

---

## 🚀 Projeyi Çalıştırma Adımları

Scriptleri aşağıda belirtilen sıra ile terminal veya RStudio üzerinden çalıştırmanız gerekmektedir:

### 1. Kütüphanelerin Kurulması
Projede kullanılan tüm paketleri (limma, WGCNA, clusterProfiler, survival, randomForest, glmnet vb.) kurmak için terminalden şu komutu çalıştırın:
```bash
Rscript packages.R
```

### 2. Verinin İndirilmesi ve Ön İşleme
Ham verileri GEO'dan indirmek, normalizasyon kontrol grafiklerini çizmek ve verileri temizleyip `data/` klasörüne kaydetmek için:
```bash
Rscript GSE68465_preprocessing.R
```
*Not: Bu komut çalıştığında eğer ham veriler yerelde yoksa NCBI GEO veritabanından indirilecektir (yaklaşık 80 MB). İnternet hızınıza bağlı olarak birkaç dakika sürebilir.*

### 3. Diferansiyel Gen İfade (DGE) Analizi
Gruplar arasında anlamlı derecede ifadesi değişen genleri limma ile tespit edip `results/GSE68465-dge.csv` olarak kaydetmek için:
```bash
Rscript GSE68465_DGE.R
```

### 4. Gen İfade Görselleştirmeleri
Genlerin dağılımlarını ve istatistiksel test sonuçlarını (Wilcoxon rank-sum test p-değerleri ile) çizdirmek için:
```bash
Rscript GSE68465_visualization.R
```
*Üretilen dosyalar:* `plots/volcanoGSE68465.png`, `plots/boxplotTop5.png`, `plots/violinTop5.png` vb.

### 5. Survival (Hayatta Kalma) Analizleri
Top 5 aday gen için Kaplan-Meier sağkalım eğrilerini çizdirmek ve klinik değişkenlerle Cox modellemesi yapmak için:
```bash
Rscript GSE68465_survival.R
```
*Üretilen dosyalar:* `plots/survival_[GEN].png` ve `plots/forest_multivariate.png`.

### 6. WGCNA Eş-ifade Ağı Analizi
Ağ analizini gerçekleştirmek, modülleri bulmak, kME değerlerine göre en iyi 5 hub geni tespit etmek ve Cytoscape ağ verilerini kaydetmek için:
```bash
Rscript GSE68465_wgcna.R
```
*Üretilen dosyalar:* `results/hubsInEachModule.csv`, `plots/module-trait2.png`, `results/edges_[modül].txt` vb.

### 7. GO Zenginleştirme (Enrichment) Analizi
Histolojik tümör derecesiyle en yüksek korelasyona sahip modüldeki genlerin işlevlerini (GO Biyolojik Süreçler) belirlemek için:
```bash
Rscript GSE68465_GO.R
```
*Üretilen dosyalar:* `plots/GO_barplot.png`, `plots/GO_dotplot.png` ve `results/GO_enrichment_results_[modül].csv`.

### 8. Makine Öğrenmesi ile Sınıflandırma (Aşama 2)
Modelleri eğitmek, veri sızıntısını önlemek için özellik seçimini fold içlerinde yapmak ve 5-katlı çapraz doğrulama ile performans skorlarını hesaplamak için:
```bash
Rscript GSE68465_ML.R
```
*Üretilen dosyalar:* `results/ML_performance_comparison.csv`, `plots/ROC_comparison.png`, `plots/confusion_matrices_all.png`.

---

## 📌 Önemli Bilgiler ve Analiz Kararları

1.  **Sınıflandırma Değişkeni (Histologic Grade):** Klinik meta-veride yer alan histolojik derece (grade) parametresi, tümörün diferansiyasyon seviyesini ("Well Differentiated" - İyi Diferansiye / "Poorly Differentiated" - Kötü Diferansiye) temsil etmektedir. Tümörlerin biyolojik agresifliğini sınıflandırmak amacıyla bu değişken ana hedef olarak seçilmiş, grade değeri bulunmayan (NA) hastalar analiz dışı bırakılmıştır.
2.  **Veri Sızıntısının Önlenmesi (Data Leakage):** Makine öğrenmesi Setup-2 (LASSO) ve Setup-3 (PCA) aşamalarında, özellik seçimi ve boyut indirgeme işlemleri çapraz doğrulama fold'larının **sadece eğitim (train) seti** üzerinde uygulanmış, test setine ait hiçbir bilgi eğitim aşamasına sızdırılmamıştır.
