# Introduction と Summary and discussion の文献比較

2026-09-14 更新。対象は Class 5・Class 6 の量子格子 Lax 演算子と、その Landau--Lifshitz 連続極限を扱う原稿。

関連する一次文献の導入・結びと、今回の局所 time-Lax 連続極限に近い構成を確認した。以下は網羅的な優先性調査ではなく、原稿で何を既知事項として認め、何を本稿の課題として残すべきかを整理するための記録である。

| 文献と確認箇所 | 既存結果 | 原稿への反映 |
| --- | --- | --- |
| M. Kruczenski, *Spin chains and string theory*, [hep-th/0311203](https://arxiv.org/abs/hep-th/0311203) | Heisenberg spin chain の coherent-state continuum limit と fast-moving string の対応。 | スピン鎖と連続 Landau--Lifshitz 理論の物理背景。 |
| J. Avan, A. Doikou, K. Sfetsos, *Systematic classical continuum limits of integrable spin chains and emerging novel dualities*, [1005.4605](https://arxiv.org/abs/1005.4605), §§2.2--2.4 | coherent-state expectation value から classical Lax/monodromy を系統的に構成。非線形演算と期待値が一般に交換しないこと、global lattice sums で coincident/overlapping indices が独立な格子和の減少により抑制されることも議論。 | coherent-state continuum Lax construction 自体を新規としない。「先行研究が overlap を無視した」とは書かない。本稿は global sum ではなく local time equation の次数評価を扱う。 |
| A. Doikou, I. Findlay, *The quantum auxiliary linear problem & Darboux--Backlund transformations*, [1706.06052](https://arxiv.org/abs/1706.06052) | closed/open lattice model の quantum time-Lax hierarchy を構成。coherent states による time evolution も議論。 | quantum time-Lax operator の一般構成自体を新規としない。finite-lattice local time equation を fixed-spin long-wavelength limit へ直接移す問題に限定する。 |
| A. Doikou, N. Karaiskos, *Generalized Landau--Lifshitz models on the interval*, [1105.5042](https://arxiv.org/abs/1105.5042) | classical continuum Hamiltonian、Lax pair、境界条件。 | continuum time matrix の classical construction を既存基盤として引用。 |
| T. Kameyama, K. Yoshida, *Anisotropic Landau--Lifshitz sigma models from q-deformed AdS₅ × S⁵ superstrings*, [1405.4467](https://arxiv.org/abs/1405.4467), §4.1 | null-like warped `SL(2)` Landau--Lifshitz model とその Lax pair。Jordanian twist と関係する非局所構造。 | draft PR #13 の独立計算により、Class 5 の traceless pair が complex field redefinition と local auxiliary gauge transformation の下でこの既知 pair に写ることを確認。新しい Class 5 classical hierarchy の発見とは扱わない。周期境界・monodromy の同一視はしない。 |
| A. Levin, M. Olshanetsky, A. Zotov, *Classical integrable systems and soliton equations related to eleven-vertex R-matrix*, [1406.2995](https://arxiv.org/abs/1406.2995) | eleven-vertex `R`-matrix に対応する classical integrable field theory / Landau--Lifshitz structure。 | Class 6 hierarchy 自体が既知であることを明示。 |
| K. Atalikov, A. Zotov, *Field theory generalizations of two-body Calogero--Moser models in the form of Landau--Lifshitz equations*, [2010.14297](https://arxiv.org/abs/2010.14297), §4 | rational eleven-vertex Landau--Lifshitz の explicit `U,V` を与える。 | 本稿 Class 6 の `U_6,V_6` と coupling/spectral/time dictionary を成分ごとに照合。これは独立検算であり新規 hierarchy の主張ではない。 |
| Ž. Krajnik, E. Ilievski, T. Prosen, V. Pasquier, *Anisotropic Landau--Lifshitz model in discrete space-time*, [2104.13863](https://arxiv.org/abs/2104.13863), Appendix A・§2 | quantum algebra の semiclassical limit から classical Sklyanin Lax matrix を得て、classical discrete-space-time zero-curvature map を構成。 | quantum-to-classical Lax limiting procedure 一般が既知であることの追加例。ただし本稿の quantum time operator の lower-symbol shared-site correction とは対象が異なる。 |
| M. de Leeuw, A. Fontanella, J. M. Nieto García, *An integrable deformed Landau--Lifshitz model with particle production?*, [2506.13598v2](https://arxiv.org/abs/2506.13598v2) | Class 5/6 continuum model、higher charges。Class 5 では quantum Lax から spatial matrix を得るが、その変数で compatible companion/time matrix は構成していない。 | Class 5 time component をその quantum lattice operator から直接導く具体的な target。pair の存在自体を新規とはしない。 |

## 現在の構成判断

Introduction では、次の三点を先に既知事項として認める。

1. coherent-state continuum Lax construction は既知。
2. quantum time-Lax hierarchy は既知。
3. Class 5/6 の最終 classical Lax structure も、局所的には既知 pair に対応する。

その上で課題を

```text
finite-lattice local quantum time equation
        -> fixed-spin long-wavelength limit
        -> surviving shared-site operator-product correction
```

に限定する。

共通設定では、第一 spatial coefficient `X` と

```text
Xi = (X^2)^downarrow - (X^downarrow)^2
```

を用い、Sutherland residual `B` を残した局所式

```text
C = 2 i (d_x Xi + [Xi,U]) - i <X1 B>
```

を提示する。対象模型では `B=0` である。Class 5、Class 6 はこの共通式の二つの具体化と、独立な classical construction への validation として置く。

Class 5 の詳細な Kameyama--Yoshida 対応、量子 L 作用素の規格化監査、実条件・周期境界条件の差は研究ノート PR #13 に保存する。paper には local `U,V` dictionary と boundary limitation のみを圧縮して載せる。

Class 6 は known eleven-vertex pair との explicit dictionary を本文に短く載せる。

XXZ は一般式の第三 control として研究ノート・検証コードに残し、第三の長い模型節にはしない。

## 現在の新規性評価

避けるべき主張は次である。

- coherent-state lower symbol の非乗法性を新しい効果と呼ぶ。
- quantum time-Lax operator の一般構成を新しい方法と呼ぶ。
- Class 5 または Class 6 の classical Lax hierarchy 自体を新しく発見したと述べる。
- Avan--Doikou--Sfetsos が overlapping-site contribution を無視したと述べる。

現時点で本稿固有の結果として残るのは、fixed physical spin `1/2` の長波長極限で **local quantum time-Lax equation** を直接展開し、二つの共有サイト積の `epsilon^2` 項の相殺後に spatial Taylor term と cubic term が `epsilon^3` に残ること、その寄与を Sutherland residual を含む共通局所恒等式として整理することである。

Class 5/6 の二模型、undeformed XXX、XXZ control、および既知 classical pair との二つの独立照合は、この局所恒等式と量子--古典 bridge の検証として位置づける。

対象文献を確認した範囲では同一の局所 shared-site formula は見つけていないが、これは literature absence や priority の証明ではない。
