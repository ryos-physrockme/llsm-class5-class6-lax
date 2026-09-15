# 共有サイト補正の解析的導出と研究ノートの追記

作業日: 2026-09-13。参照コミット: `722b5bf14bb64ca895375ab05569e11dd69f9712`。

## 今回確かめたこと

有限格子の時間Lax方程式の共有サイト積について、coherent-state lower symbolの差を解析的に導いた。第一空間係数をX、そのlower symbolをUとし、第二モーメントの差をXi=(X^2)↓−U^2と定義する。局所物理空間は二次元、状態はrank-one射影の積、Hamiltonianの先頭項は交換演算子Pを用いた2Pとする。

三重項部分空間で一次Hamiltonian補正がスカラーになる条件の下で、共有サイト補正は

```text
C = 2i(∂x Xi + [Xi,U]) − i <X1 B>
```

となる。BはSutherland relationの格子間隔二次係数の残差、< >は同一点の二サイト積状態のsymbolである。B=0なら局所時間行列の補正Delta V=2i Xi+f(lambda)Iに吸収できる。f(lambda)は位置とスピンに依存しないスカラーである。第二空間係数を含む各項の相殺と、任意の一次Hamiltonian補正で残る追加項も明示した。

さらに、Uのスピン3成分への依存が階数3なら、三重項の条件は独立に仮定する必要がない。Sutherland残差のサイト交換に対して偶な部分から、一次Hamiltonian補正の偶成分が全スピンと可換であることが従い、三重項・一重項の重複度が一であることから条件を導ける。Class 5/6では全ての非零スペクトルパラメータで必要な小行列式が非零となる。階数条件を外すとこの推論が使えない例も検証した。

導出は補助空間の行列成分を非可換なまま行っており、補助空間を二次元に限らない。物理空間が二次元である条件は維持する。古典Hamiltonian流との同定、全運動方程式と曲率の同値性、有限時間の量子発展の誤差評価は別問題である。

## 研究ノートと検証コード

- `notes/parts/10_shared_site_generalization.tex`: 成立条件、二次項のTaylor展開、三次項の四つの寄与、Sutherland残差、階数条件による仮定の削減、模型の適用範囲を記述。
- `notes/full_research_note_ja.tex`: 上記節を索引の直前に組み込み。既存本文の数式は変更していない。
- `notes/shared_site_generalization_ja.tex`: 同一節から独立した抜粋PDFを作るためのソース。
- `scripts/generalization/verify_shared_site_proof.py`: 非可換な補助空間成分を保持する独立実装。
- `results/generalization/shared_site_proof.json`: 今回の34チェックの実行結果とソースハッシュ。

```bash
make verify-generalization
make note-generalization
make note
```

既存の50チェックと今回の34チェックをローカルで実行した。計84項目のうち80項目が厳密にゼロ、4項目が仮定を外した検算で予期した非零となった。この数は独立な研究結果の数ではなくassertion数である。環境はPython 3.13.5、SymPy 1.14.0。

今回の独立実装のSHA-256:
`24571fb4c63c31626900c8daa6646a354793643b5a604031d5de64edf4e3062e`。

既存の全模型検証と機械学習の再実行はローカルでは行っていない。既存CIが呼ぶ`verify-class6`に`verify-generalization`を依存させ、PRで追加検証も実行するようにした。抜粋PDFはLuaLaTeXでローカル作成し、式番号、引用、ページ配置を確認した。完全ノートのPDFはPRのCIビルドで確認する。作業ブランチには生成PDFをコミットしない。

## 文献と評価

Avan–Doikou–Sfetsos, arXiv:1005.4605の第2節と、Doikou–Findlay, arXiv:1706.06052の第2・5節を比較対象として確認した。前者は期待値と非線形演算が一般に交換しない点も扱う。後者は量子時間Lax成分の階層を構成する。これらの枠組み自体を新規とは主張しない。本節の局所的な格子間隔三次補正式について、文献上の優先性は未確定である。

確認した一次資料に基づき、ノートの参考文献の3種類の書誌情報を訂正した。Corcoranの頭文字をL.に、arXiv:1406.2995の題名と頁を *Classical integrable systems and soliton equations related to eleven-vertex R-matrix*, Nucl. Phys. B 887 (2014) 400–422に、arXiv:1706.06052の題名と掲載情報を *The quantum auxiliary linear problem & Darboux-Backlund transformations*, PoS CORFU2019 (2020) 210に直した。互換性のため既存の引用キーは維持する。他の参考文献全体を再点検したという意味ではない。

一次資料: https://arxiv.org/abs/1005.4605 、https://arxiv.org/abs/1706.06052v3 、https://arxiv.org/abs/1406.2995 、https://arxiv.org/abs/2306.10423 。PDFのスクリーンショット取得はツール内部エラーとなったため、文献比較は抽出テキストと書誌ページに基づく。

論文原稿`paper/`は変更していない。一般式の新規性調査と、Class 6と既知eleven-vertex古典Lax対の直接照合は未完了であり、今回の成功のみで出版上の評価を確定しない。
