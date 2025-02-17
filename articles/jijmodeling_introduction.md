---
title: "代数的数理モデラー JijModelingの紹介"
emoji: "👋"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["数理最適化"]
published: false
---

こんにちは、株式会社Jijの山城です。

この記事は[Jij Inc. Advent Calendar 2024](https://adventar.org/calendars/10963)の2日目の記事です。

Jijでは数理最適化向けのソフトウェアの開発を行っています。その中でも、数理最適化問題をモデリングするためのライブラリ JijModeling について紹介します。

# JijModelingとは

JijModelingは数理最適化問題をモデリングするためのライブラリです。数理最適化問題を解くためには、問題を数学的にモデル化する必要があります。JijModelingは、そのような数理最適化問題をモデル化するための機能を提供します。

JijModelingを使うことで、数理最適化問題を簡単にモデル化することができます。例えば、以下のようなコードで、最小化問題をモデル化することができます。

```python
import jijmodeling as jm

# 変数の定義
x = jm.

