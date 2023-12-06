---
title: "Simulated Quantum Annealingを導出する"
emoji: "🦁"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["physics", "quantum", "annealing"]
published: false
publication_name: jij_inc
---


初めまして、株式会社Jijの山城です。

この記事はJij Inc. Advent Calendar 2023の7日目の記事です。  
前回はJijの広報による ["ChatGPTを使って広報活動レポートを作成し業務効率化"](https://zenn.dev/jij_inc/articles/371d0736797b3e) でした。

Jij(ジェイアイジェイ)は数理最適化、量子技術を用いたソフトウェアの開発を行っています。


Jijの社名の由来となっているのは統計物理のイジングモデルという数理モデルです。Jijはイジングモデルの基底状態を求める量子アルゴリズムである量子アニーリングを技術シードとして設立されたスタートアップです。現在では量子アニーリングだけでなくより広く数理最適化技術によるソフトウェア開発を行っていますが、7日目のアドベントカレンダーではJijらしく量子アニーリングに関することを書きます。

量子アニーリングはもちろん量子デバイスで実行されるアルゴリズムなのですが、よくシミュレーションする手段として Simulated Quantum Annealing (SQA) という手法が用いられます。この記事ではSQAを導出することでSQAがどのような手法なのかを紹介します。

# はじめに

## 想定読者

物理学科の学部生くらいで、基礎的な量子力学、統計力学を学んでいることを想定します。簡単なブラケット演算や分配関数などは既知とします。

## 量子アニーリング

量子アニーリングを模したアルゴリズムとしてシミュレーテッド量子アニーリングというアルゴリズムがあります。量子アニーリングを通常のコンピュータでシミュレーションする方法の一つです。量子アニーリングは量子力学つまり、シュレディンガー方程式に従って実行されるのですが、シュレディンガー方程式を愚直にシミュレートするのは通常のコンピュータでは困難です。なのでシュレディンガー方程式を使わず平衡統計力学の手法を使ってシミュレートします。

しかし平衡統計力学とシュレディンガー方程式によるダイナミクスはそもそも異なる物理です。しかしなぜ平衡統計力学を持ち出すのでしょうか？大雑把には量子アニーリングにおいて断熱時間発展を考えて時間発展中もハミルトニアンの基底状態付近にいるとすると、低温の平衡状態（基底状態）で近似してもよさげな結果になるだろうという考えから来ています。

この記事では横磁場のみの量子アニーリングに対するシミュレーテッド量子アニーリング (Simulated Quantum Annealing; SQA) を紹介します。特にアルゴリズムを導出するのに必要なトロッター展開を用いた量子モンテカルロ法を中心とした記事です。

SQAを実際に実行してみたい方はOpenJijをぜひ利用してみてください。

# SQAの導出

## Hamiltonian of the transverse field quantum annealing

横磁場量子アニーリングのハミルトニアンは以下で記述されます。

$$\hat H (t) = A(t) \sum_i \hat \sigma_i^x + B(t) \left(\sum_i \hat \sigma_i^z + \sum_{i,j}J_{ij}\hat \sigma^z_i \hat \sigma^z_j\right) = A(t) \sum_i \hat \sigma_i^x + B(t)H_p(\hat \sigma^z)$$
where
$$\sum_i \hat \sigma_i^z + \sum_{i,j}J_{ij}\hat \sigma^z_i \hat \sigma^z_j = H_p (\hat \sigma^z)$$
このハミルトニアン下での分配関数を考えましょう。

## Partition function

平衡統計力学でのこのハミルトニアン下での状態を考えるわけなので分配関数を考えましょう。

$$Z(\beta) = \mathrm{Tr}\exp\left( -\beta \hat H(t) \right) = \mathrm{Tr} \exp\left(A(t) \sum_i \hat \sigma_i^x + B(t)H_p(\hat \sigma^z)\right)$$

実は以下の計算の手続きでこの$\mathrm{Tr}$ を処理して古典的なスピン配位の和に変換することで、量子系の分配関数を古典系の分配関数へ書き換えることができます。

そうすることで古典系の分配関数になってしまえばあとはそれをそのままモンテカルロ法すればよいわけです。シミュレーテッドアニーリングがモンテカルロ法をベースであることを思い出すと、この系に対してのモンテカルロ法をそのままアニーリングのようにパラメータを変化させながら実行すれば量子アニーリングっぽいアルゴリズムが作れることが想像できると思います。では早速この分配関数に対しての計算を進めていきましょう。


## Suzuki-Trotter decomposition

そのままでは指数の中身の2つのハミルトニアンが非可換なので計算を進めることができません。ここで鈴木トロッター展開と呼ばれる方法を使って指数関数を2つに分解します。

$$Z =\mathrm{Tr} \exp\left(A(t) \sum_i \hat \sigma_i^x + B(t)H_p(\hat \sigma^z)\right)\\
= \mathrm{Tr} \left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right)\exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right]^m + \mathcal O\left(\left(\frac{\beta}{m}\right)^2\right)$$

あるスピン配位を $\vec\sigma$ と書くことにします。トレースはこのスピン配位の全パターンでの状態和なので、

$$Z \simeq \sum_{\vec \sigma} 
\bra{\vec \sigma}
\left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right)\exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right]^m \ket{\vec \sigma}$$

です。中の演算は同じ形の行列の指数関数が$m$個ありますが、その積の間に完全系の式

$$\hat 1 = \hat 1_k = \sum_{\vec \sigma} | \vec \sigma_k\rang \lang \vec \sigma_k|$$


を挟みます。ここで添字$k$は$\hat 1$をたくさん入れていくので区別がつくように入れました。分配関数に$\hat 1$を挟んでいくと

$$
Z \simeq \sum_{\vec \sigma}
\left \lang\vec \sigma\left|
\hat 1_1
\left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right) \exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right] \hat 1_2 \left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right)\exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right]\hat 1_3\cdots
\right| \vec \sigma\right\rang
$$

$$
= \sum_{\vec \sigma}
\left \lang \vec \sigma\left| \sum_{\vec \sigma_1}| 
\vec \sigma_1\rang \lang 
\vec \sigma_1|
\left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right) \exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right] \sum_{\vec \sigma_2}| \vec \sigma_2\rang \lang \vec \sigma_2|\cdots
\right| \vec \sigma\right\rang\\
=\prod_{k=1}^m \sum_{\vec \sigma_k}\left \lang \vec \sigma_k\left| 
\left[\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right) \exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right)\right]
\right| \vec \sigma_{k+1}\right\rang
$$

となります。ここで $m+1 = 1$ とします。また

$$
\exp\left(-\frac{\beta}{m}B(t)H_p(\hat \sigma^z)\right) |\vec \sigma_k\rang = \exp\left(-\frac{\beta}{m}B(t)H_p(\vec \sigma_k)\right) |\vec \sigma_k\rang
$$

なので$H_p$に関するexpは行列からただの数に変わります。よってそのまま外にくくりだします。次に

$$
\lang\vec \sigma_k|\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right) |\vec \sigma_{k+1}\rang
$$

を考えます。あとはこの横磁場項の一体問題を解けばよいです。スピン配位を1スピンごとのテンソル積 $\ket{\vec \sigma_k} = \otimes_i \ket{(\sigma_i)_k}$ で記述すると

$$
\lang\vec \sigma_k|\exp\left(-\frac{\beta}{m}A(t) \sum_i \hat \sigma_i^x\right) |
\vec \sigma_{k+1}\rang = \prod_{i} \lang(\sigma_i)_k|\exp\left(-\frac{\beta}{m}A(t) \hat \sigma_i^x\right) |(\sigma_i)_{k+1}\rang
$$

となります。ではこの一体問題を考えましょう。行列の指数関数はテイラー展開で定義されるのでそのままexpを展開してやって、$(\hat \sigma_i^x)^2 = 1$であることを使うと、

$$
\exp\left(-\frac{\beta}{m}A(t) \hat \sigma_i^x\right) = \hat 1 \cosh\left(-\frac{\beta A(t)}{m}\right) + \hat \sigma_i^x \sinh\left(-\frac{\beta A(t)}{m}\right)
$$

となります。よって

$$
\lang(\sigma_i)_k|\exp\left(-\frac{\beta}{m}A(t) \hat \sigma_i^x\right) |(\sigma_i)_{k+1}\rang = \lang(\sigma_i)_k|\left(\hat 1 \cosh\left(-\frac{\beta A(t)}{m}\right) + \hat \sigma_i^x \sinh\left(-\frac{\beta A(t)}{m}\right)\right) |(\sigma_i)_{k+1}\rang \\
= \cosh\left(-\frac{\beta A(t)}{m}\right) \lang(\sigma_i)_k|(\sigma_i)_{k+1}\rang +  \sinh\left(-\frac{\beta A(t)}{m}\right) \lang(\sigma_i)_k|-(\sigma_i)_{k+1}\rang
$$

となります。ここで$\hat \sigma_i^x$はスピンを反転させるので$\hat \sigma_i^x|\sigma_i\rang = |-\sigma_i\rang$を使いました。

ここまでで分配関数は

$$
Z \simeq \prod_k\sum_{\vec \sigma_k}\exp\left(-\frac{\beta}{m}B(t)H_p(\vec \sigma_k)\right) \left(\cosh\left(-\frac{\beta A(t)}{m}\right) \lang(\sigma_i)_k|(\sigma_i)_{k+1}\rang +  \sinh\left(-\frac{\beta A(t)}{m}\right) \lang(\sigma_i)_k|-(\sigma_i)_{k+1}\rang\right)
$$

となりました。$\hat \sigma^x$ からの寄与も指数の肩に載せて整理します。

ここで、

$$
\cosh\left(-\frac{\beta A(t)}{m}\right) = \sqrt{\cosh\left(-\frac{\beta A(t)}{m}\right) \sinh\left(-\frac{\beta A(t)}{m}\right)\left(\tanh\left(-\frac{\beta A(t)}{m}\right)\right)^{-1}}
=\sqrt{\frac{1}{2}\sinh\left(-2\frac{\beta A(t)}{m}\right)\left(\tanh\left(-\frac{\beta A(t)}{m}\right)\right)^{-1}}\\
\sinh\left(-\frac{\beta A(t)}{m}\right) = \sqrt{\frac{1}{2}\sinh\left(-2\frac{\beta A(t)}{m}\right)\tanh\left(-\frac{\beta A(t)}{m}\right)}
$$


であることを使います。この時 $\cosh$ と $\sinh$ の違いは $\tanh$ のべきだけなので以下のように $\hat \sigma^x$ からの寄与を計算することができます。整理すると

$$
\cosh\left(-\frac{\beta A(t)}{m}\right)\lang(\sigma_i)_k|(\sigma_i)_{k+1}\rang +  \sinh\left(-\frac{\beta A(t)}{m}\right) \lang(\sigma_i)_k|-(\sigma_i)_{k+1}\rang\\
= \sqrt{\frac{1}{2}\sinh\left(-2\frac{\beta A(t)}{m}\right)\tanh\left(-\frac{\beta A(t)}{m}\right)^{-(\sigma_i)_k (\sigma_i)_{k+1}}}\\
=\sqrt{\frac{1}{2}\sinh\left(-2\frac{\beta A(t)}{m}\right)}\exp\left(-\frac{1}{2}\ln\tanh\left(-\frac{\beta A(t)}{m}\right){(\sigma_i)_k (\sigma_i)_{k+1}}\right)
$$

となります。よって分配関数は

$$
Z \simeq \prod_k\sum_{\vec \sigma_k}
\exp\left(-\frac{\beta}{m}\left(
B(t)H_p(\vec \sigma_k)
-\frac{m}{2\beta}\ln\tanh\left(-\frac{\beta A(t)}{m}\right){(\sigma_i)_k (\sigma_i)_{k+1}}
\right)\right)\\
=\sum_{\vec \sigma_1,\cdots,\vec \sigma_m}
\exp\left(-\frac{\beta}{m}\left(
B(t)\sum_k H_p(\vec \sigma_k)
-\frac{m}{2\beta}\ln\tanh\left(-\frac{\beta A(t)}{m}\right)\sum_k{(\sigma_i)_k (\sigma_i)_{k+1}}
\right)\right)
$$

となります。スピン配位にかかわらない係数は省略しました。これで帰着できたこの指数関数の中身はもとのスピン変数に一つ添え字が増えて $k, k+1$ 間に相互作用が追加された系となっています。

# 最後に
＼Rustエンジニア・数理最適化エンジニア募集中！／
株式会社Jijでは、数学や物理学のバックグラウンドを活かし、量子計算と数理最適化のフロンティアで活躍するRustエンジニア、数理最適化エンジニアを募集しています！
詳細は下記のリンクからご覧ください。皆さんのご応募をお待ちしております！
Rustエンジニア: https://open.talentio.com/r/1/c/j-ij.com/pages/51062
数理最適化エンジニア: https://open.talentio.com/r/1/c/j-ij.com/pages/75132

