# 局所 time-Lax 連続極限と共有サイト補正の文献比較

調査日: 2026-09-13。
対象ブランチ: `research/shared-site-generalization`。

本メモは、`notes/parts/10_shared_site_generalization.tex` で導出した局所補正式の文献上の位置づけと、Class 5 / Class 6 の既知 Lax 構造との関係を確認するための調査記録である。既存の coherent-state symbol calculus、量子 time-Lax 演算子、連続 Lax 構成そのものを新規と扱わない。ここでいう「確認できなかった」は、指定した文献と検索範囲で同一式を見つけなかったという意味であり、文献上の不在を証明するものではない。

## 1. Avan--Doikou--Sfetsos (2010)

J. Avan, A. Doikou and K. Sfetsos,
*Systematic classical continuum limits of integrable spin chains and emerging novel dualities*,
Nucl. Phys. B 840 (2010) 469--490,
[arXiv:1005.4605](https://arxiv.org/abs/1005.4605).

この論文は量子スピン鎖から古典連続 Lax 行列を構成する一般手順を与えている。したがって、coherent-state expectation value を用いた量子鎖から連続 Lax 理論への移行自体は既存である。

特に第2節では、積状態における monodromy の期待値がサイトごとの Lax 演算子の期待値へ因数分解されることを使い、

\[
L_{a i}=1+\delta\,l_{a i}+O(\delta^2)
\]

から連続 monodromy を構成する。展開中のある項が \(\delta^n\) を伴い、独立なサイト和が \(M\) 個あれば、その項の continuum scaling は独立和の個数も含めて判定される。この power counting が global な monodromy / conserved Hamiltonians の議論を支えている。

同論文第2.4節では、一般に

\[
\langle F(A)\rangle\neq F(\langle A\rangle)
\]

であることを明記している。また、高次の local Hamiltonians の多項式について、局所因子のサイト添字が coincident または overlap する項は generic な非重複項より独立なサイト和が一つ少なくなるため、global continuum limit では相対的に抑制されると論じている。したがって、「既存研究が expectation value の非乗法性や overlapping sites を無視した」という説明は不正確である。

今回の計算との違いは対象にある。今回展開しているのは global な格子和ではなく、各サイトで成り立つ局所零曲率式

\[
\partial_{t_{\rm lat}}L_n=A_{n+1}L_n-L_nA_n
\]

である。ここには、独立なサイト和の個数が減ることによる抑制機構はない。規格化後に

\[
\widehat L^\downarrow=1+\epsilon U+O(\epsilon^2),\qquad
\widehat A^\downarrow=\epsilon^2 V^\downarrow+O(\epsilon^3),
\qquad t=\epsilon^2t_{\rm lat}
\]

とすると、共有サイト積の個々の connected contribution は \(O(\epsilon^2)\) から始まる。二つの向きの同一点値は相殺するが、隣接状態の最初の Taylor 係数と同一点の三次係数が \(O(\epsilon^3)\) に残る。この次数は continuum time evolution と同じである。

したがって、本稿で強調できる差は

- global conserved quantities における overlap の suppression ではなく、
- local quantum time-Lax equation における overlap の surviving contribution

である。

ADS (2010) を対立する先行研究として扱うべきではなく、空間 Lax の continuum construction と power counting の既存基盤として引用した上で、局所 time equation では異なる次数評価が必要になることを説明するのが妥当である。

## 2. Doikou--Findlay (2017/2020)

A. Doikou and I. Findlay,
*The quantum auxiliary linear problem & Darboux--Backlund transformations*,
PoS CORFU2019 (2020) 210,
[arXiv:1706.06052](https://arxiv.org/abs/1706.06052).

この論文は closed / open integrable lattice models に対して、量子 Lax pair の time components の hierarchy を量子 monodromy から系統的に構成している。したがって、有限格子の量子 time-Lax 演算子の一般的構成そのものは今回の新規事項ではない。

第5節では local operator の量子時間発展を扱い、第5.2節で coherent states と path integral を導入して semi-classical description を議論する。ただし、今回確認した範囲では、量子 time-Lax 演算子を fixed local spin の長波長極限へ取り、局所零曲率式中の \(A_{n+1}L_n\) と \(L_nA_n\) の共有物理サイトに由来する lower-symbol correction を continuum time order まで展開する計算は行っていない。

したがって、量子 time-Lax hierarchy と coherent-state dynamics は既知の二つの材料である一方、今回の研究対象はそれらを

\[
\text{finite-lattice local time equation}
\longrightarrow
\text{fixed-spin long-wavelength continuum time matrix}
\]

として直接つなぐ際の局所積の扱いに限定すべきである。

## 3. de Leeuw--Fontanella--Nieto García (2026): Class 5

M. de Leeuw, A. Fontanella and J. M. Nieto García,
*An integrable deformed Landau--Lifshitz model with particle production?*,
SciPost Phys. 20 (2026) 041,
[arXiv:2506.13598](https://arxiv.org/abs/2506.13598).

Class 5 の coherent-state continuum model、Hamiltonian flow、保存量、および量子 Lax 演算子から得られる classical spatial Lax matrix はこの論文で既に構成されている。したがって、Class 5 continuum model や spatial Lax matrix の存在は今回の新規事項ではない。

一方、同論文は spatial Lax matrix の Poisson algebra を確認した後、運動方程式を抽出する compatible companion matrix を見つけられなかったと明記している。したがって、今回の Class 5 計算は「新しい模型を発見した」という位置づけではなく、この未完成だった time component を有限格子の time operator から導出し、さらに独立な classical r-matrix / monodromy construction と照合したものと位置づけられる。

## 4. eleven-vertex Landau--Lifshitz hierarchy と Class 6 の直接照合

A. Levin, M. Olshanetsky and A. Zotov,
*Classical integrable systems and soliton equations related to eleven-vertex R-matrix*,
Nucl. Phys. B 887 (2014) 400--422,
[arXiv:1406.2995](https://arxiv.org/abs/1406.2995),
および
K. Atalikov and A. Zotov,
*Field theory generalizations of two-body Calogero--Moser models in the form of Landau--Lifshitz equations*,
J. Geom. Phys. 164 (2021) 104161,
[arXiv:2010.14297](https://arxiv.org/abs/2010.14297)
を比較した。

Class 6 の量子模型が eleven-vertex model に対応すること、および eleven-vertex classical r-matrix が rational Landau--Lifshitz hierarchy を与えることは既知である。Atalikov--Zotov (2021), Sec. 4 は、その rational model の空間・時間 Lax 行列を標準基底で明示しているため、今回の Class 6 の最終 \(U_6,V_6\) と成分ごとに比較できる。

Class 6 のスピンから

\[
\Sigma=\begin{pmatrix}S^3&S^+\\S^-&-S^3\end{pmatrix},
\qquad \Sigma^2=1
\]

を定義する。非零複素数 \(\beta\) を

\[
\beta^2=-2\alpha_6
\]

で定義し、既知模型のスペクトルパラメータと空間微分係数を

\[
z=\beta\lambda,\qquad k=-\frac4\beta
\]

と置く。Class 6 空間行列の traceless part は厳密に

\[
U_{6,0}(\lambda)
=\frac{\beta}{4}U_{\rm 11v}(\beta\lambda,\Sigma)^{\mathsf T}
\]

となる。Class 6 の完全な \(U_6\) との差は \(1/(4\lambda)\) 倍の単位行列であり、零曲率を変えない。

時間規約を

\[
\tau=\frac{i\beta^2}{2}t=-i\alpha_6 t
\]

と変換すると、time matrix も

\[
V_6(\lambda)
=\frac{i\beta^2}{2}V_{\rm 11v}(\beta\lambda,\Sigma)^{\mathsf T}
=-i\alpha_6 V_{\rm 11v}(\beta\lambda,\Sigma)^{\mathsf T}
\]

として全成分が一致する。転置は補助空間の基底規約の差に対応し、交換子の順序を反転する。係数 \((\beta/4)k=-1\) により空間微分項を含む零曲率規約も一致する。

さらに rational eleven-vertex model の anisotropy map を用いると、同じ \(\beta,z,k,\tau\) の辞書で運動方程式も Class 6 の三成分式へ一致する。この一致は `scripts/class6/verify_eleven_vertex_comparison.py` で、スピン成分と一階・二階空間微分を独立な記号として検証する。

この結果は Class 6 の hierarchy の新規性を支持するものではない。むしろ、有限格子の量子 time-Lax 演算子から本研究で得た pair が、既知の classical eleven-vertex hierarchy に正しく着地するという独立検算である。Class 5 における independent classical r-matrix / monodromy comparison と同じ役割を Class 6 にも与える。

## 5. 現段階での新規性評価

今回の条件付き一般式

\[
\mathcal C
=2i\bigl(\partial_x\Xi+[\Xi,U]\bigr)
-i\langle X_1\mathcal B\rangle,
\qquad
\Xi=(X^2)^\downarrow-U^2
\]

と、Sutherland relation の該当次数 \(\mathcal B=0\) における

\[
\mathcal C=2i\bigl(\partial_x\Xi+[\Xi,U]\bigr)
\]

について、上記の一次文献では同一の局所 time-Lax continuum formula を確認できなかった。ただし、これは targeted literature audit の結果であり、優先性の証明ではない。

少なくとも次の主張は避ける必要がある。

1. coherent-state lower symbol の非乗法性を新しい効果と呼ぶこと。
2. quantum time-Lax operator の一般的構成を新しい方法と呼ぶこと。
3. ADS (2010) が shared-site / overlapping-index contributions を無視したと述べること。
4. Class 6 の Landau--Lifshitz hierarchy 自体を新しいと述べること。

一方、原稿の中心候補として残るのは次である。

- fixed spin \(1/2\) の長波長極限において、有限格子の局所 time-Lax equation を直接 continuum へ移す際の operator-product correction の power counting を明示すること。
- その補正が Sutherland relation の残差を含む局所恒等式として整理でき、指定した条件下で adjoint covariant derivative 型の局所 time correction へ吸収できること。
- Class 5 では既存研究で未構成だった time component を得ること。
- Class 5 では classical r-matrix / monodromy construction、Class 6 では既知 rational eleven-vertex Landau--Lifshitz pair と、量子格子から得た continuum time matrix がそれぞれ独立に一致すること。
- XXZ を既知模型の第三 control として、一般式が二つの非Hermitian模型だけの偶然ではないことを確認すること。

ここまで進むと、論文の構成は単に「Class 5 と Class 6 を計算した」より、

\[
\text{local quantum time equation}
\longrightarrow
\text{conditional continuum identity}
\longrightarrow
\text{Class 5 / Class 6 applications and independent validations}
\]

と置く方が、実際の内容を正確に表す。ただし、この一般恒等式自体が標準的な symbol calculus と Sutherland relation から比較的短く導けることも事実であり、概念的新規性を過大評価しない。

## 6. 原稿への反映方針

一般式を原稿に入れる場合、Introduction では方法の新規性を広く主張せず、既存の二つの流れを先に認めるべきである。

1. coherent-state continuum Lax construction: ADS (2010)。
2. quantum time-Lax hierarchy: Doikou--Findlay (2017/2020) など。

その上で、「本稿では finite-lattice local time equation を fixed-spin long-wavelength limit へ直接移し、その局所積が continuum time order に残る条件と形を調べる」と課題を限定する。

Section 2 には、Sutherland relation の残差 \(\mathcal B\)、第二モーメント差 \(\Xi\)、一般式を導入し、解析的証明の主要部分を示す候補がある。Class 5 / Class 6 の各節は、その一般式の具体化と独立な Hamiltonian / Lax check として再配置できる。Class 5 の classical r-matrix comparison と Class 6 の eleven-vertex comparison を対称的な validation として見せる価値がある。XXZ は第三の長い模型節ではなく、一般式の既知模型による短い control とする。

この段階では paper source をまだ変更せず、研究ノートと検証コードに結果を固定する。次の原稿改訂では、一般式を Section 2 にどの長さで入れるか、二模型の重複した導出をどこまで圧縮できるかを先に設計し、既存の詳細検算を失わない構成にする。
