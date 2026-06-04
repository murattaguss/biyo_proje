# Biyoenformatik Proje Raporu: Akciğer Adenokarsinomu (GSE68465) Analizi

## 1. Projenin Amacı ve Özeti
Bu projenin amacı, Akciğer Adenokarsinomu (Lung Adenocarcinoma) hastalarına ait gen ifade verilerini biyoinformatik yaklaşımlarla analiz etmektir. Çalışmada, tümör dokularının histolojik diferansiyasyon dereceleri (Histologic Grade) baz alınmıştır. Bu kapsamda, iyi diferansiye (Well Differentiated) tümörler ile kötü diferansiye (Poorly Differentiated) tümörler moleküler düzeyde karşılaştırılmıştır.

Proje temel olarak iki ana aşamadan oluşmaktadır. Birinci aşamada, NCBI GEO veritabanından elde edilen Affymetrix ham mikrodizin verileri RMA (Robust Multi-array Average) algoritmasıyla normalize edilmiştir. İki farklı grup (Well vs Poorly) arasında ifadesi değişen genleri bulmak için Diferansiyel Gen İfadesi (DGE) analizi uygulanmış, seçilen genler üzerinden Kaplan-Meier ve Cox regresyon sağkalım analizleri yapılmıştır. Ağırlıklı Gen Eş-İfade Ağı Analizi (WGCNA) ile tümör derecesiyle ilişkili gen modülleri ve hub genler belirlenmiş, ardından Gen Ontolojisi (GO) analiziyle bu modüllerin işlevleri incelenmiştir. İkinci aşamada ise farklı özellik seçimi yöntemleri (biyolojik, LASSO, PCA) kullanılarak Random Forest, Destek Vektör Makineleri (SVM) ve Lojistik Regresyon algoritmalarıyla makine öğrenmesi modelleri kurulmuş ve hastalığın diferansiyasyon seviyesini tahmin etme performansları karşılaştırılmıştır.

---

## 2. Veriseti Seçimi ve Ön İşleme

### 2.1. Veri Seti Seçimi ve Klinik Özellikler
Çalışmada, NCBI GEO veritabanından alınan ve akciğer adenokarsinomu (LUAD) vakalarını içeren GSE68465 veri seti kullanılmıştır. Veri setinin klinik yapısı incelendiğinde 443 tümör ve 19 normal doku örneği içerdiği belirlenmiştir. Kanser evresinin doğrudan sınıflandırma etiketi olarak seçilmesi kısıtlandığından, analizlerde hedef parametre olarak histolojik derece (histologic grade) seçilmiştir.

Veri ön işleme adımlarında grade bilgisi bulunmayan normal örnekler (19 adet) ile geçersiz veriye sahip tümör örnekleri (7 adet) çıkarılmış ve 436 tümör örneği üzerinden analize devam edilmiştir. Kalan hastaların histolojik derece dağılımı şu şekildedir:
- İyi Diferansiye (Well Differentiated): 60 hasta
- Orta Diferansiye (Moderate Differentiation): 209 hasta
- Kötü Diferansiye (Poorly Differentiated): 167 hasta

Veri setinde ayrıca hayatta kalma (survival) analizleri için gerekli olan yaş, cinsiyet, sigara geçmişi, son temas veya ölüm süresi (months to last contact or death) ve yaşam durumu (vital status: Alive/Dead) verileri eksiksiz bulunmaktadır.

### 2.2. Veri Normalizasyonu ve Gen Anotasyonu
Ham `.CEL` formatındaki dosyalar R ortamına alınarak RMA (Robust Multi-array Average) yöntemiyle normalize edilmiştir. Bu adımda arka plan gürültüsü düzeltilip veriler log2 tabanına dönüştürülmüştür. Yapılan işlemin doğruluğunu göstermek için rastgele seçilen 30 hastanın normalizasyon öncesi ve sonrası dağılımları boxplot grafiği ile gösterilmiştir.

* **Figür 1:** Normalizasyon öncesi ve sonrası verilerin boxplot karşılaştırması (`plots/normalization_before_after.png`).

Normalizasyon işleminden sonra, mikrodizin çipindeki prob kodları `hgu133a.db` paketi kullanılarak geçerli gen sembollerine (Gene Symbol) çevrilmiştir. Karşılığı bulunmayan veya boş olan problar çıkarılmıştır. Birden fazla probun tek bir gene denk geldiği durumlarda, varyansı en yüksek olan prob genin temsilcisi olarak korunup diğerleri silinmiştir. Matrisin satırları ve sütunları yer değiştirilerek veriseti Hasta × Gen formatına getirilmiştir. Hesaplama sürelerini kısaltmak ve gürültüyü azaltmak için verideki varyansı en yüksek olan 3000 gen (Top-3000 Highly Variable Genes) filtrelenerek DGE dahil tüm sonraki analiz adımlarına aktarılmıştır.

---

## 3. Diferansiyel Gen İfadesi (DGE) Analizi
İki uç fenotipin karşılaştırması için DGE analizine İyi Diferansiye (Well, n=60) ve Kötü Diferansiye (Poorly, n=167) grupları dahil edilmiş, Orta Diferansiye (Moderate, n=209) grup dışarıda bırakılmıştır. Bu iki grup arasındaki moleküler farklılıkları belirlemek amacıyla `limma` paketiyle lineer modelleme yürütülmüştür. Tasarım matrisinde (design matrix) kesim noktası (intercept) sıfırlanarak doğrudan grupların kendi gen ifade katsayıları modellenmiştir. İstatistiksel anlamlılık sınırı olarak kat değişim eşiği |LogFC| > 1 ve çoklu hipotez testi düzeltmeli yanlış keşif oranı FDR (Adjusted P Value) < 0.05 kullanılmıştır. Bayesyen eBayes düzeltmesi sonucunda elde edilen diferansiyel gen ifadesi özeti Tablo 1'de gösterilmiştir.

| Karşılaştırma Grubu | Down-Regulated (Azalan İfade) | Not Significant (Anlamsız) | Up-regulated (Artan İfade) | Toplam Anlamlı DEG |
| :--- | :---: | :---: | :---: | :---: |
| **Poorly vs. Well Differentiated** | 89 Gen | 2833 Gen | 78 Gen | 167 Gen |

* **Tablo 1:** Well vs. Poorly Diferansiye LUAD Tümörleri Arasındaki DGE Dağılımı.

Analiz sonucunda saptanan 167 anlamlı Diferansiyel İfade Edilen Gen (DEG), tümörün kötü diferansiye faza geçerken sergilediği transkriptomik değişimleri yansıtmaktadır. Toplam 78 gen Poorly grubunda yukarı regüle (Up) olurken, 89 gen ise aşağı regüle (Down) olmuştur. Bu dağılımın gen genelindeki istatistiksel anlamlılık (adj.P.Val) ve kat değişimi (logFC) ilişkisi Volcano Plot grafiği ile özetlenmiştir.

* **Figür 2:** Diferansiyel ifade gösteren genlerin dağılımını gösteren Volcano Plot (`plots/volcanoGSE68465.png`).

Kat değişimi (logFC) mutlak değerine göre sıralanan genler arasında en yüksek farka sahip ilk 5 gen (Top-5 DEG) şunlardır: **`PGC`**, **`CYP4B1`**, **`C1orf116`**, **`ADH1B`** ve **`SFTPC`**. Bu genlerin tamamının Poorly (Kötü diferansiye) grupta azaldığı (down-regulated) görülmüştür. Bu genler içinden en belirgin değişimi gösteren gen `PGC` genidir (logFC = -2.65, adj.P.Val = 9.06 $\times 10^{-20}$). `PGC` geninin iki grup arasındaki ekspresyon dağılımı, Wilcoxon Rank-Sum test sonucu eklenmiş Boxplot ve Violin grafikleriyle görselleştirilmiştir.

* **Figür 3:** `PGC` genine ait histolojik derece bazlı boxplot dağılımı (`plots/boxplot_PGC.png`).
* **Figür 4:** `PGC` genine ait histolojik derece bazlı violin dağılımı (`plots/violin_PGC.png`).

Analizde en yüksek değişim gösteren ilk 5 genin (`PGC`, `CYP4B1`, `C1orf116`, `ADH1B`, `SFTPC`) tamamını içeren toplu Violin ve Boxplot şemaları da `facet_wrap` fonksiyonu ile oluşturulmuştur.

* **Figür 5:** Top-5 DEG genini içeren toplu boxplot grafikleri (`plots/boxplotTop5.png`).
* **Figür 6:** Top-5 DEG genini içeren toplu violin grafikleri (`plots/violinTop5.png`).

---

## 4. Sağkalım Analizi ve Cox Regresyon Modeli
DGE analizinde keşfedilen genlerin hastaların genel sağkalımıyla ilişkili olup olmadığını test etmek amacıyla sağkalım analizi uygulanmıştır. Proje yönergesine uygun olarak, DGE listesinden dengeli bir set oluşturulmuştur: 3 adet down-regüle gen (`PGC`, `CYP4B1`, `C1orf116`) ve 2 adet up-regüle gen (`COL11A1`, `RRM2`). Hastalar, her bir genin ekspresyon medyan değerine göre "High" (Yüksek İfade) ve "Low" (Düşük İfade) olmak üzere iki gruba ayrılmış ve Kaplan-Meier sağkalım eğrileri çizilmiştir.

* **Figür 7:** 5 aday gen için oluşturulan Kaplan-Meier sağkalım eğrileri (`plots/survival_PGC.png`, `plots/survival_CYP4B1.png`, `plots/survival_C1orf116.png`, `plots/survival_COL11A1.png`, `plots/survival_RRM2.png`).

Gen ifadelerinin yanı sıra hastanın yaşı, cinsiyeti, sigara geçmişi ve tümör derecesi gibi klinik değişkenlerin de sağkalım üzerindeki etkisini değerlendirmek amacıyla **Çok Değişkenli Cox Orantılı Hazard Regresyon Modeli** (Multivariate Cox Proportional Hazards Model) kurulmuştur. Model, eksiksiz klinik veriye sahip 435 vaka üzerinden eğitilmiştir.

Kurulan çok değişkenli Cox modelinin genel uyum kabiliyeti (goodness of fit) Concordance Index (C-Index) = 0.637 (se = 0.02) olarak ölçülmüştür. Likelihood Ratio Test sonucu $p = 2 \times 10^{-5}$ gelerek modelin istatistiksel açıdan anlamlı olduğu görülmüştür. Model parametrelerinin incelenmesinden elde edilen başlıca bulgular aşağıdadır:

* **Yaş (Age):** Hazard Ratio (HR) = 1.0310, %95 Güven Aralığı (CI) [1.0168 - 1.0450] ve p-değeri = $1.55 \times 10^{-5}$ (***). Diğer tüm değişkenler kontrol altındayken, yaşın her bir yıl artışının ölüm riskini bağımsız olarak %3.1 oranında artırdığı bulunmuştur.
* **Cinsiyet (Sex - Male):** HR = 1.3701, %95 CI [1.0384 - 1.8080] ve p-değeri = 0.0260 (*). Erkek hastaların ölüm tehlikesinin (hazard) kadın hastalara kıyasla %37.0 daha yüksek olduğu tespit edilmiştir.
* **C1orf116 (Yüksek İfade Grubu):** HR = 0.7390, %95 CI [0.5367 - 1.0180] ve p-değeri = 0.0638 (.). Bu genin yüksek grupta yer almasının, ölüm tehlikesini %26.1 oranında azaltıcı bir eğilim gösterdiği gözlemlenmiştir.
* **RRM2 (Yüksek İfade Grubu):** HR = 1.3645, %95 CI [0.9947 - 1.8720] ve p-değeri = 0.0540 (.). Bu genin tümör dokusundaki yüksek ifadesinin, ölüm riskini bağımsız olarak %36.4 oranında artırma eğiliminde olduğu görülmüştür.

Bu regresyon analizinin katsayıları ve güven aralıkları Forest Plot ile özetlenmiştir.

* **Figür 8:** Çok değişkenli Cox regresyon modeline ait katsayıları ve Hazard Ratio (HR) değerlerini gösteren Forest Plot (`plots/forest_multivariate.png`).

---

## 5. Ağırlıklı Gen Eş-İfade Ağı Analizi (WGCNA)
Top-3000 Highly Variable Gen (HVG) matrisi kullanılarak, genlerin birbirleriyle olan fonksiyonel ilişkilerini incelemek amacıyla Ağırlıklı Gen Eş-İfade Ağı Analizi (WGCNA) yürütülmüştür. Analizin ilk adımında "Ölçekten Bağımsız Topoloji" (Scale-Free Topology) modeline uygunluğu sağlamak için `pickSoftThreshold` fonksiyonu ile en uygun güç parametresi aranmıştır. Model uyum ($R^2$) değerinin 0.8 eşiğini aştığı ilk güç parametresi olan **softPower = 6** seçilmiştir.

* **Figür 9:** WGCNA analizi için ölçekten bağımsız topoloji modeli uyum ve ortalama bağlantısallık grafikleri (`plots/soft_th.png`).

Seçilen softPower=6 değeri kullanılarak Topological Overlap Matrix (TOM) oluşturulmuş ve hiyerarşik kümeleme algoritması ile gen dendrogramı elde edilmiştir. Dinamik dal kesimi (dynamic tree cut) ile genler, modül başına en az 30 gen düşecek şekilde renkli modüllere atanmıştır. Korelasyonu yüksek (cutHeight = 0.25) modüller birleştirilerek toplam 6 modül (blue, brown, grey, red, turquoise, yellow) tanımlanmıştır. Standart WGCNA pratiğine uygun olarak hiçbir modüle atanamayan grey modül, işlevsel analizlerin dışında tutulmuştur.

* **Figür 10:** Hiyerarşik kümeleme ağacı ve dinamik dal kesimiyle belirlenen gen modülleri (`plots/gene_dendrogram_module.png`).
* **Figür 11:** Modül eigengene'lerinin kümelenmesi ve birleştirilme sınırı grafiği (`plots/METree_merged.png`).

Modüllerin klinik niteliklerle (yaş, sigara, hayatta kalma ve histolojik derece) ilişkisini göstermek için Pearson korelasyon matrisi hesaplanarak bir Module-Trait İlişkisi Heatmap şeması oluşturulmuştur.

* **Figür 12:** Gen modülleri ile klinik parametreler arasındaki Pearson korelasyonlarını gösteren Module-Trait Relationships Heatmap (`plots/module-trait2.png`).

Analiz sonucunda hastanın ölüm durumu (`vital_status`, Pearson $r = 0.23$) ve kötü tümör derecesi (`grade_poorly`, Pearson $r = 0.40$) ile en güçlü ilişkiyi gösteren modülün **Turquoise (Turkuaz)** modülü olduğu saptanmıştır. Bu durum, modülün hastalık ilerlemesiyle pozitif ilişki gösteren modül olarak değerlendirilmesini sağlamıştır.

Proje kuralı gereğince, her modülün en yüksek bağlantısallığa sahip Top-5 çoklu hub geni belirlenerek kaydedilmiştir. Turkuaz modülüne ait ağ yapısı, 0.05 TOM eşiğinde filtrelenerek `edges_turquoise.txt` ve `nodes_turquoise.txt` formatlarında dışa aktarılmıştır.

* **Figür 13:** Turkuaz modülündeki genlerin circular layout düzeniyle oluşturulan gen etkileşim ağı ve hub genlerin (sarı dikdörtgenler) gösterimi (`plots/cytofull.png`). Ağ görseli harici Cytoscape uygulamasında oluşturulmuştur.

Turkuaz modülünün en yüksek üyeliğe (kME) sahip en önemli 5 hub geni şunlardır: **`MAD2L1`** (kME = 0.929), **`CEP55`** (kME = 0.900), **`BUB1B`** (kME = 0.900), **`RACGAP1`** (kME = 0.896) ve **`PRC1`** (kME = 0.895).

---

## 6. Gen Ontolojisi (GO) Zenginleştirme Analizi
WGCNA analizinde hastalık fenotipiyle en yüksek korelasyonu gösterdiği saptanan Turquoise (Turkuaz) modülündeki genlerin, ilişkili olduğu hücresel mekanizmalar ve biyolojik süreçleri incelemek amacıyla Gen Ontolojisi (GO) analizi yapılmıştır. Gen sembolleri doğrudan SYMBOL formatında `enrichGO()` fonksiyonuna aktarılmıştır. Analizde biyolojik süreçler (Biological Process - BP) kategorisi kullanılmış, çoklu test düzeltmesi için Benjamini-Hochberg (BH) yöntemi tercih edilmiştir. Anlamlılık sınırları p-value < 0.05 ve q-value < 0.2 olarak belirlenmiştir.

* **Figür 14:** Turkuaz modül genlerine ait en anlamlı GO Biyolojik Süreçlerini gösteren Barplot (`plots/GO_barplot.png`).
* **Figür 15:** Turkuaz modül genlerine ait en anlamlı GO Biyolojik Süreçlerini gösteren Dotplot (`plots/GO_dotplot.png`).

### Biyolojik Süreçlerin Fonksiyonel Yorumu
GO zenginleştirme analizinin sonuçları incelendiğinde, Turkuaz modülündeki genlerin ağırlıklı olarak **Kardeş Kromatid Ayrılması** (Sister Chromatid Segregation, p = $1.16 \times 10^{-23}$), **Mitotik Çekirdek Bölünmesi** (Mitotic Nuclear Division, p = $6.65 \times 10^{-22}$), **Kromozom Ayrılması** (Chromosome Segregation, p = $1.34 \times 10^{-23}$), **DNA Replikasyonu** (DNA Replication, p = $1.96 \times 10^{-18}$) ve **Hücre Döngüsü Faz Geçişinin Düzenlenmesi** (Regulation of Cell Cycle Phase Transition, p = $9.80 \times 10^{-16}$) süreçlerinde görev aldığı saptanmıştır. 

Bu durum, akciğer tümör dokusu kötü diferansiye (Poorly) faza geçerken hücresel çoğalma (proliferasyon) ve bölünme mekanizmalarıyla ilişkili biyolojik süreçlerin aktifleştiğini göstermektedir. Modüldeki `RRM2` ve `BUB1B` gibi hub genlerin hücre döngüsündeki işlevleri, tümörün agresif bölünme eğilimiyle ve düşük sağkalımla uyumlu bir tablo ortaya koymaktadır.

---

## 7. Makine Öğrenmesi Sonuçları

### 7.1. Özellik Seçimi ve Çıkarımı Stratejileri (Setup Düzenekleri)
Çalışmanın ikinci aşamasında, hastaların transkriptomik profilleri kullanılarak tümörün histolojik derecesini (Well vs. Poorly) tahmin edebilen sınıflandırma modelleri geliştirilmiştir. Farklı özellik yönetimi yaklaşımlarını karşılaştırmak amacıyla 3 farklı "Setup" düzeneği kurulmuştur:

* **Setup-1 (Biyolojik Özellik Seti):** DGE analizinden elde edilen en anlamlı 50 DEG gen ile WGCNA analizinden elde edilen modüllerin lider hub genleri birleştirilerek rafine bir biyolojik özellik seti oluşturulmuştur.
* **Setup-2 (Algoritmik Özellik Seçimi - LASSO):** Biyolojik ön bilgi kullanılmaksızın, en değişken 3000 gen `glmnet::cv.glmnet()` fonksiyonu üzerinden L1-regularizasyonuna (LASSO) sokulmuştur. Veri sızıntısını (data leakage) önlemek amacıyla bu işlem yalnızca eğitim (Train) verisi üzerinde çapraz doğrulamayla uygulanarak katsayısı sıfır olmayan özellikler seçilmiştir.
* **Setup-3 (Algoritmik Özellik Çıkarımı - PCA):** Verideki varyansı özetlemek amacıyla Temel Bileşen Analizi (PCA) uygulanmıştır. Yine veri sızıntısını engellemek için dönüşüm matrisi yalnızca eğitim setinde hesaplanmış, ardından en güçlü ilk 10 Temel Bileşen (PC1-PC10) özellik olarak kullanılmıştır.

### 7.2. Model Performanslarının Karşılaştırılması ve Değerlendirme Metrikleri
Oluşturulan 3 farklı Setup düzeneği, üç farklı makine öğrenmesi algoritması olan Random Forest (RF), Destek Vektör Makineleri (SVM) ve Lojistik Regresyon (GLM) modellerinde test edilmiştir. Modellerin hiperparametreleri ve performansı Stratified 5-Fold Cross-Validation altında değerlendirilerek sonuçlar raporlanmıştır (Tablo 2).

| Deney No | Özellik Seti (Setup) | Sınıflandırma Modeli | Accuracy | Precision | Recall | F1 Score | AUC |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| 1 | Setup-1 (Biyolojik) | Random Forest (RF) | 0.9031 | 0.9240 | 0.9461 | 0.9349 | 0.9597 |
| 2 | **Setup-1 (Biyolojik)** | **SVM** | **0.9163** | **0.9405** | **0.9461** | **0.9433** | **0.9460** |
| 3 | Setup-1 (Biyolojik) | Lojistik Regresyon (GLM) | 0.8370 | 0.9167 | 0.8563 | 0.8854 | 0.8621 |
| 4 | **Setup-2 (LASSO)** | **Random Forest (RF)** | **0.9119** | **0.9249** | **0.9581** | **0.9412** | **0.9635** |
| 5 | **Setup-2 (LASSO)** | **SVM** | **0.9163** | **0.9353** | **0.9521** | **0.9436** | **0.9493** |
| 6 | Setup-2 (LASSO) | Lojistik Regresyon (GLM) | 0.8855 | 0.9379 | 0.9042 | 0.9207 | 0.9051 |
| 7 | Setup-3 (PCA) | Random Forest (RF) | 0.8811 | 0.9023 | 0.9401 | 0.9208 | 0.9368 |
| 8 | Setup-3 (PCA) | SVM | 0.9075 | 0.9345 | 0.9401 | 0.9373 | 0.9505 |
| 9 | Setup-3 (PCA) | Lojistik Regresyon (GLM) | 0.8943 | 0.9333 | 0.9222 | 0.9277 | 0.9376 |

* **Tablo 2:** Sınıflandırma Modellerinin ve Setup Düzeneklerinin Performans Karşılaştırma Matrisi.

### 7.3. ML Sonuçlarının Yorumlanması ve Değerlendirilmesi
Elde edilen sonuçlara göre en yüksek sınıflandırma doğruluğuna **%91.63 Accuracy** değeri ile **Setup-1 (Biyolojik Özellikler) ve Setup-2 (LASSO) üzerinde eğitilen SVM modelleri** ulaşmıştır. Modellerin ROC eğrisi performansı baz alındığında ise **0.9635 AUC** ve **%91.19 Accuracy** oranı ile **Setup-2 (LASSO) düzeneğinde eğitilen Random Forest (RF) modeli** dikkat çekmektedir.

* **Biyolojik Özelliklerin Gücü (Setup-1 vs. Setup-2):** Biyolojik tabanlı özellik seti (Setup-1) ile LASSO algoritmik özellik seçim (Setup-2) düzeneklerinin birbirine çok yakın sonuçlar vermesi, DGE ve WGCNA yöntemleriyle elde edilen genlerin tümör derecesini yansıtmada veri güdümlü yöntemlerle kıyaslanabilir bir ayırt edicilik sergilediğini ortaya koymaktadır.
* **PCA Boyut İndirgeme (Setup-3):** İlk 10 temel bileşenin kullanıldığı Setup-3 düzeneğinde veri boyutu düşürülmüş olup, buna karşın SVM modelinde %90.75 doğruluk oranlarına ulaşılarak PCA'in etkili bir boyut indirgeme alternatifi olduğu görülmüştür.

Sınıflandırıcıların ROC eğrileri ve 3x3 grid formatındaki hata matrisleri (confusion matrices) analiz edilerek kaydedilmiştir.

* **Figür 16:** 3 farklı setup düzeneği için modellerin ROC eğrisi karşılaştırmaları (`plots/ROC_comparison.png`).
* **Figür 17:** Tüm model ve setup kombinasyonlarına ait 3x3 grid hata matrisleri (`plots/confusion_matrices_all.png`).

---

## 8. Kaynakça

1. Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47.
2. Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559.
3. Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284-287.