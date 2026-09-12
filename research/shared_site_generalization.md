# スピン1/2鎖のcoherent-state連続極限における時間Lax行列の補正

計算日: 2026-09-13。原稿の参照版: `ac97ec22c21c5f0c8feb5ca3aaf73e1a72ea2325`。

本ノートは、最近接相互作用を持つ量子スピン鎖の時間Lax演算子にcoherent-state lower symbolを適用する際の局所補正を検討する。既存原稿の二模型に共通する式について、成立の十分条件と、その条件を外した場合の残差を求める。以下の一般恒等式とXXZ鎖への適用は今回の計算であり、文献から引用した結果としては扱わない。文献上の優先性は未確認である。原稿本文と既存の検証プログラムは変更していない。

## 1. 設定と十分条件

物理的な局所Hilbert空間と補助空間をそれぞれ二次元とする。二つの物理サイトを1、2で表し、その交換演算子をPとする。補助空間の添字は省略し、X_1、X_2などの添字は演算子が作用する物理サイトを指定する。

格子間隔を$\epsilon$、量子スペクトルパラメータを$u$、連続理論のスペクトルパラメータを$\lambda$とし、$u=\lambda/\epsilon$を用いる。規格化した空間Lax演算子と二サイトHamiltonianを

\[
\widehat L=1+\epsilon X+\epsilon^2Y+O(\epsilon^3),\qquad
h=2P+\epsilon h^{(1)}+O(\epsilon^2),
\tag{1}
\]

と展開する。X、Yは補助空間と一つの物理サイトに作用する4×4行列であり、原稿の第一、第二空間係数に相当する。さらに

\[
\partial_u\widehat L=\epsilon^2Z+O(\epsilon^3)
\tag{2}
\]

でZを定義する。今回のスペクトルスケーリングでは$Z=\partial_\lambda X$であり、微分時には連続結合定数を固定する。一般恒等式の検証ではX、Y、Zをまず独立な行列として扱う。

第一の条件は、Sutherland relation

\[
[h,\widehat L_2\widehat L_1]
=\widehat L_2\partial_u\widehat L_1
 -(\partial_u\widehat L_2)\widehat L_1
\tag{3}
\]

の$\epsilon$二次係数が成立することである。その左辺から右辺を引いた残差を

\[
\mathcal B:=
[2P,X_2X_1+Y_2+Y_1]
+[h^{(1)},X_2+X_1]-Z_1+Z_2
\tag{4}
\]

と定義する。条件は$\mathcal B=0$である。Pは物理サイトを交換するため、[P,Y_1+Y_2]=0は任意のYについて成立する。

第二の条件は、三重項部分空間への射影$P_{\rm sym}=(1+P)/2$を用いて

\[
P_{\rm sym}h^{(1)}P_{\rm sym}=cP_{\rm sym}
\tag{5}
\]

となる、スピンに依存しないスカラーcが存在することである。この条件により、同じ向きのスピンの積状態でのh^(1)の期待値はcとなる。スピンの向きに依存するepsilon一次のポテンシャルがepsilon二次の交換エネルギーを支配しない、という今回の連続極限と整合する十分条件である。必要条件としては主張しない。

二次元の物理空間では、この条件を満たすh^(1)全体は

\[
h^{(1)}=h_-+c\,1+d(1-P),\qquad Ph_-P=-h_-
\tag{6}
\]

と書ける。c、dは任意の複素スカラーである。三重項部分空間に直交する一重項部分空間が一次元であることを用いた。

第三の条件は、各サイトでrank-one射影

\[
\rho(x)=\frac{1+S^i(x)\sigma^i}{2},\qquad
\rho^2=\rho,\quad\operatorname{tr}\rho=1,\quad S^iS^i=1
\tag{7}
\]

を用いることである。$\sigma^i$はPauli行列であり、$i=1,2,3$について和を取る。スピンの積は複素双線形に拡張してよい。隣接サイトには滑らかなrho(x)を割り当てる。これは大スピン極限や有限時間の量子発展の近似を仮定するものではない。

## 2. 演算子積から得られる恒等式

一サイト演算子のlower symbolを$X^\downarrow=\operatorname{tr}_{\rm phys}(\rho X)$、二サイト演算子の同一点でのsymbolを$\langle Q\rangle=\operatorname{tr}_{12}[(\rho\otimes\rho)Q]$と定義する。補助空間の行列積の順序は維持する。連続空間行列Uと第二モーメントの差Xiを

\[
U=X^\downarrow,\qquad \Xi=(X^2)^\downarrow-U^2
\tag{8}
\]

と定義する。

規格化した時間Lax演算子は

\[
\widehat A=-i\widehat L_2^{-1}
\bigl([h,\widehat L_2]+\partial_u\widehat L_2\bigr)
=\epsilon A^{(1)}+\epsilon^2A^{(2)}+O(\epsilon^3),
\tag{9}
\]

であり、その二係数は

\[
\begin{aligned}
A^{(1)}&=-i[2P,X_2],\\
A^{(2)}&=iX_2[2P,X_2]-i[2P,Y_2]
-i[h^{(1)},X_2]-iZ_2.
\end{aligned}
\tag{10}
\]

原稿と同じく、重なった積のsymbolから個別symbolの積を引いた量を

\[
\begin{aligned}
D^+_n&=(\widehat A_{n+1}\widehat L_n)^\downarrow
       -\widehat A_{n+1}^\downarrow\widehat L_n^\downarrow,\\
D^-_n&=(\widehat L_n\widehat A_n)^\downarrow
       -\widehat L_n^\downarrow\widehat A_n^\downarrow
\end{aligned}
\tag{11}
\]

と定義する。D^+の物理サイトはn,n+1、D^-の物理サイトはn-1,nである。epsilon展開の二次、三次係数をD_2^±、D_3^±と書く。

任意のXについて、同一点で

\[
D_2^+(\rho,\rho)=D_2^-(\rho,\rho)=2i\Xi
\tag{12}
\]

となる。したがって二次項は差D^+-D^-で相殺する。$\rho_x=\partial_x\rho$を用いた最初の空間Taylor係数の和は

\[
\delta_{\rho_R}D_2^+(\rho,\rho)[\rho_x]
+\delta_{\rho_L}D_2^-(\rho,\rho)[\rho_x]
=2i\partial_x\Xi .
\tag{13}
\]

ここでdeltaは指定した左または右の行列引数に関する方向微分である。二つが加算されるのは、D^-自体を引くこととrho(x-epsilon)の展開の負号による。式(12)、(13)にはSutherland relationもh^(1)の条件も不要である。

三次係数には条件が関わる。式(5)のもとで、**Sutherland relationの残差をまだゼロとせずに**、次の恒等式を得た。

\[
D_3^+(\rho,\rho)-D_3^-(\rho,\rho)
=2i[\Xi,U]-i\langle X_1\mathcal B\rangle .
\tag{14}
\]

したがって、$\epsilon^{-3}(D^+-D^-)$の極限を$\mathcal C$とすれば

\[
\mathcal C=2i\bigl(\partial_x\Xi+[\Xi,U]\bigr)
-i\langle X_1\mathcal B\rangle .
\tag{15}
\]

式(4)のSutherland relationが成立する場合には、最後の項が消える。結果として原稿の共通式

\[
\mathcal C=\partial_x\Delta V+[\Delta V,U],\qquad
\Delta V=2i\Xi+f(\lambda)1
\tag{16}
\]

を回収する。f(lambda)は位置とスピンに依存しない任意のスカラーである。連続時間$t$と格子時間$t_{\rm lat}$の関係$t=\epsilon^2t_{\rm lat}$のもとで、この$\mathcal C$が時間発展と同じ次数に現れる。

### 恒等式の確認方法と証明の範囲

式(10)を式(11)に代入する。二次項はA^(1)X、三次項はA^(2)XおよびA^(1)Yから生じる。交換演算子の恒等式$PX_1=X_2P$、および$P(\rho\otimes\rho)=(\rho\otimes\rho)P=\rho\otimes\rho$を使うと、式(12)が従う。$\rho$の微分には$\rho\rho_x+\rho_x\rho=\rho_x$を用いる。

式(14)は同じ演算子積から直接検証している。任意のrank-one射影を定数の物理基底変換で$\operatorname{diag}(1,0)$に移し、その一般の接ベクトルを$\left(\begin{smallmatrix}0&p\\q&0\end{smallmatrix}\right)$と置く。ここで$p,q$は独立な複素数である。X、Y、Zにはそれぞれ独立な16個の記号を置き、h^(1)には式(6)の一般形を置く。このとき式(14)の差の全成分が恒等的にゼロになる。これはスピン方向やLax係数の数値サンプリングではなく、全成分が未定の行列に関する多項式恒等式の検証である。定数の同時物理基底変換はPと式(5)を保存するので、任意のrank-one射影とその接ベクトルに戻せる。

三次項中の$Y$依存は、$A^{(2)}$に含まれる$-i[2P,Y_2]$と$A^{(1)}Y$の寄与の間で相殺する。したがって一般式はYを独立な追加データとして必要としない。特定模型で(X^2)にYが現れる恒等式は、同じ情報の別の表示である。この点は原稿の「第二空間係数が寄与する」という説明と矛盾しないが、第二係数を常に独立の物理情報として強調するのは避けるべきである。

実際、$X=X_0\otimes1+X_i\otimes\sigma^i$と展開し、各$X_i$を補助空間行列と定義すると、$U=X_0+X_iS^i$である。Pauli行列の積から

\[
\Xi=\sum_{i,j=1}^3
\frac{\partial U}{\partial S^i}
\frac{\partial U}{\partial S^j}
\left(\delta_{ij}-S^iS^j+i\varepsilon_{ijk}S^k\right)
\tag{16a}
\]

とも表せる。ここで$\varepsilon_{123}=1$で、$k$についても和を取る。行列積の順序は表示通りである。したがって、第一係数のsymbolである$U$のスピン依存性全体が既知なら、補正に必要な第二モーメントの差も復元できる。これはPauli代数の帰結であり、新しいsymbol calculusの主張ではない。

さらに、同一点での直接時間symbolには

\[
(A^{(2)})^\downarrow+2i\Xi
=-i\langle[h^{(1)},X_2]\rangle-iZ^\downarrow
\tag{17}
\]

が成り立つ。Xのsymbolがスピンにaffineであることを使えば、補正後の時間行列は

\[
V=-2(\boldsymbol S\times\partial_x\boldsymbol S)^i
       \frac{\partial U}{\partial S^i}
-i\langle[h^{(1)},X_2]\rangle-iZ^\downarrow+f(\lambda)1
\tag{18}
\]

と書ける。式(18)は有限時間の量子発展の収束や、任意の模型での曲率と運動方程式の同値性まで証明するものではない。これらは別の問題である。

## 3. 既存の二模型とXXZ鎖での検算

Class 5とClass 6は、de Leeuwらのスピン1/2鎖の分類に由来する二つの非Hermitian変形である。今回の検証では原稿の演算子を具体的に

\[
K_5=\begin{pmatrix}0&1&-1&0\\0&0&0&0\\0&0&0&0\\0&0&0&0\end{pmatrix},\qquad
K_6=\begin{pmatrix}0&1&1&0\\0&0&0&-1\\0&0&0&-1\\0&0&0&0\end{pmatrix}
\tag{19}
\]

と定義する。基底の順番は$\{|\uparrow\uparrow\rangle,|\uparrow\downarrow\rangle,|\downarrow\uparrow\rangle,|\downarrow\downarrow\rangle\}$である。

Class 5では、連続結合$\alpha_5$を使って$X=P/(2\lambda)+\alpha_5 K_5$、$Y=0$、$h^{(1)}=2\alpha_5 K_5$とする。$PK_5P=-K_5$なので式(5)を満たす。Class 6では、連続結合$\alpha_6$を使って$X=P/(2\lambda)+\alpha_6\lambda K_6+\alpha_6^2\lambda^3 K_6^2$、$Y=\alpha_6 K_6/2+\alpha_6^2\lambda^2 K_6^2$、$h^{(1)}=0$とする。両模型とも$Z=\partial_\lambda X$である。両者について式(4)、(12)--(18)を独立に検算し、原稿の時間行列をスカラー規格化も含めて回収した。

第三例として既知のXXZ鎖を使う。格子異方性eta、量子スペクトルパラメータuに対し、正則R行列を

\[
R(u;\eta)=\frac1{\sinh\eta}
\begin{pmatrix}
\sinh(\eta+2\eta u)&0&0&0\\
0&\sinh(2\eta u)&\sinh\eta&0\\
0&\sinh\eta&\sinh(2\eta u)&0\\
0&0&0&\sinh(\eta+2\eta u)
\end{pmatrix}
\tag{20}
\]

と定義する。R(0)=Pで、h=P partial_u R(0)の交換係数はeta→0で2になる。これは標準的なsix-vertex R行列のスペクトル規格化を変えたものである。今回のコードでは指数変数q=exp(eta)、t=exp(2eta u)を使い、正則性と量子Yang--Baxter方程式も記号計算で確認した。

連続異方性$\gamma$を$\eta=\epsilon\gamma$で定義し、$u=\lambda/\epsilon$、$\widehat L=(\sinh\eta/\sinh(2\eta u))R$を用いる。コードとJSONでは$\gamma$を`g`と表す。スペクトル関数を$a=\gamma\coth(2\gamma\lambda)$、$b=\gamma\operatorname{csch}(2\gamma\lambda)$と定義すると

\[
X=\begin{pmatrix}a&0&0&0\\0&0&b&0\\0&b&0&0\\0&0&0&a\end{pmatrix},\qquad
Y=\frac{\gamma^2}{2}\operatorname{diag}(1,0,0,1),\qquad h^{(1)}=0.
\tag{21}
\]

$a'=-2b^2$、$b'=-2ab$、$a^2-b^2=\gamma^2$（プライムは$\lambda$微分）を用いると式(4)と補正式が成立する。スピンに依存しない定数を除いたHamiltonian密度は

\[
\mathcal H=-\frac12\partial_x\boldsymbol S\cdot\partial_x\boldsymbol S
+\frac{\gamma^2}{2}(S^3)^2.
\tag{22}
\]

Poisson bracketを$\{S^i(x),S^j(y)\}=2\varepsilon_{ijk}S^k(x)\delta(x-y)$とする。$\delta(x-y)$はDiracのデルタ関数、$\varepsilon_{ijk}$は式(16a)の完全反対称テンソルである。このとき

\[
\partial_t\boldsymbol S=-2\boldsymbol S\times
\bigl(\partial_x^2\boldsymbol S+\gamma^2S^3\boldsymbol e_3\bigr),\qquad
\boldsymbol e_3=(0,0,1)^{\mathsf T}.
\tag{23}
\]

式(18)の時間行列はこの運動方程式とoff-shellで対応することを確認した。運動方程式(23)の左辺から右辺を引いた残差の成分を$E_i$とすると、曲率は$(bE_1\sigma^1+bE_2\sigma^2+aE_3\sigma^3)/2$となり、$ab\ne0$のスペクトル領域では逆関係も取れる。これは既知XXZ模型の検算であり、新しいXXZ模型や新しい可積分階層の主張ではない。

## 4. 条件を外した検算

任意の演算子Xで式(16)が成立するわけではない。例えば

\[
X=\sigma^1\otimes\sigma^1+\sigma^3\otimes\sigma^3,\qquad
Y=Z=h^{(1)}=0,\qquad \rho=\operatorname{diag}(1,0),\quad\rho_x=0
\tag{24}
\]

では$U=\sigma^3$、$\Xi=1$であり、共通式が予言する$2i[\Xi,U]$はゼロである。しかし直接の三次項は$-4i\sigma^3$となる。その差は式(15)の$-i\langle X_1\mathcal B\rangle$で正確に説明される。これはSutherland条件を落とせない例であり、可積分鎖の中での反例ではない。

## 5. 原稿への対応の優先順位

| 対応 | 判断 | 理由 |
| --- | --- | --- |
| 既存のcontinuum-Lax構成との比較 | 必須 | 2010年の先行研究は量子と連続極限の区別や期待値の非乗法性にも言及している。「以前は単純にfactorizeしていた」とは書かず、本稿で追跡する局所時間方程式とepsilon三次項を特定する。 |
| 条件付き一般式の導出 | 有力な追加候補 | 二例の共通パターンから、仮定を指定した恒等式に進んだ。独立の式変形確認と重複文献調査を経て、共通設定か比較節に短く入れる価値がある。 |
| XXXでの係数不一致の説明 | 維持し、導入でも短く示す | 本稿の直接symbolを用いる手順では補正が必要であることの明確な検算。既存のXXX理論への新しい量子効果とは主張しない。 |
| XXZの詳細を本文に追加 | 必須ではない | 短い例または補遺で十分。新しい長い模型節を追加するより、一般式の適用例として使う。 |
| Class 6とeleven-vertex古典Lax対の直接対応表 | 追加推奨だが今回未実施 | 既知構造との対応の主張を読み手が検証できるようにする。一般恒等式とは独立した作業である。 |
| MLのAbstractからの削除 | 必須ではない | 補助的な係数探索であることが明確なら、残すかどうかは編集上の判断。計算の妥当性を左右しない。 |
| 題名変更、大スピン極限、高次Hamiltonian、有限時間誤差評価 | 今回は不要 | 題名は現在の限定を明記すれば維持可能。他は追加研究の範囲であり、この局所恒等式の条件に混ぜない。 |

これにより「Class 5/6でたまたま同じ形になっただけ」という懸念は弱まる。一方、「既存のsymbol calculusとSutherland relationから導ける恒等式が、独立の新規成果としてどの程度評価されるか」は残る。今回の成功だけで出版上の価値や新規性を確定してはいけない。

## 6. 文献上の位置づけ

以下は今回実際に確認した一次資料である。一般恒等式の優先性を証明する網羅的調査ではない。

- J. Avan, A. Doikou, K. Sfetsos, *Systematic classical continuum limits of integrable spin chains and emerging novel dualities*, [arXiv:1005.4605v2](https://arxiv.org/abs/1005.4605v2), §§2.2--2.4, 4. 式(2.14)--(2.22)はmonodromyの積状態期待値とpower countingを扱い、§2.4は期待値と非線形演算が交換しないことを明示している。式(4.6)--(4.10)はXXZのR行列と連続空間行列を扱う。今回の局所時間方程式のepsilon三次補正と同じ計算をしている、とまでは確認できていない。
- A. Doikou, I. Findlay, *The quantum auxiliary linear problem & Darboux-Backlund transformations*, [arXiv:1706.06052](https://arxiv.org/abs/1706.06052), §§2, 5. 量子時間Lax行列の階層を構成し、coherent statesによる時間発展の扱いにも触れる。本稿の有限格子時間演算子の背景として比較すべき文献であり、この枠組み全体を新規とは主張できない。
- M. de Leeuw, A. Fontanella, J. M. Nieto Garcia, *An integrable deformed Landau-Lifshitz model with particle production?*, [arXiv:2506.13598v2](https://arxiv.org/abs/2506.13598v2), §2, §4, Appendix A. Class 5/6の模型と連続スケーリングの出発点である。

PDF資料は抽出テキストを確認した。webのPDF screenshotは呼び出したが内部エラーとなったため、画像による照合は完了していない。

## 7. 再現方法と制限

実行プログラム: `scripts/generalization/verify_shared_site_identity.py`。
結果: `results/generalization/shared_site_identity.json`。

```bash
python scripts/generalization/verify_shared_site_identity.py
```

Python 3.13.5、SymPy 1.14.0で50チェックが成功した。一般行列の恒等式、Class 5/6の既存時間行列、XXZの正則性・Yang--Baxter方程式・スケーリング・off-shell関係、および条件を外した非零の残差を含む。JSONには実行環境、検証対象、スクリプトのSHA-256を記録する。

既存リポジトリの全検証スイートとMLの再学習は実行していない。主張は二次元の補助・物理空間、rank-oneの積状態symbol、式(5)とSutherland relationの指定次数に限る。別の局所スピン表現、高次Hamiltonian、有限時間の量子発展の収束は今回確認していない。
