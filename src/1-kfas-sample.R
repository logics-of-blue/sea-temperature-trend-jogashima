
# ライブラリ読み込み ---------------------------------------------------------------

library(KFAS)
library(dplyr)
library(lubridate)
library(forecast)

packageVersion("KFAS")
packageVersion("dplyr")
packageVersion("lubridate")
packageVersion("forecast")

# データ読み込み -----------------------------------------------------------------

# CSVファイルの読み込み(シミュレーションで生成された模擬データ)
sample_data <- read.csv(file = "ssm-sample-data.csv")

# time　　　：日付(論文では0.5か月だが、本サンプルコードでは月単位データとした)
# temp　　　：水温
# distance　：黒潮流軸までの距離
# kuroshio_a：黒潮流路がA型のときに1をとるフラグ
head(sample_data)

# 秋冬の黒潮Aフラグの追加
sample_data <- 
  sample_data %>% 
  mutate(winter = as.numeric(month(as.Date(time)) %in% c(11, 12, 1, 2, 3)),
         winter_a = winter * kuroshio_a)

# winter　：秋冬に1をとるフラグ
# winter_a：秋冬かつ、黒潮がAであるときに1をとるフラグ
head(sample_data)

# ts型に変換
# 0.5か月単位データの場合はfrequency = 24と指定する
sample_ts <- 
  sample_data %>% 
  select(-time, -winter) %>%             # 日付列と秋冬フラグを削除
  ts(start = c(1998, 1), frequency = 12) # ts型に変換(1998年1月開始。12か月1周期)

# 1998年12月までのデータを取得
window(sample_ts, end = c(1998, 12))

# 時系列折れ線グラフ
autoplot(sample_ts, facets = TRUE)


# 訓練データと検証データに分割 ----------------------------------------------------------

# 訓練データ
train <- window(sample_ts, end = c(2021, 12))

# 検証データ
validate <- window(sample_ts, start = c(2022, 1))


# モデル化を行う関数の作成 --------------------------------------------------------------------

# モデルを作る関数(黒潮流路のデータ利用)
make_ssm_kuroshio <- function(ts_data) {
  # モデルの構造を決める
  build_ssm <- SSModel(
    H = NA,
    temp ~
      SSMtrend(degree = 2,                  # 平滑化トレンドモデル
               Q = c(list(0), list(NA))) +
      SSMseasonal(
        sea.type = "dummy", # ダミー変数を利用した季節成分
        period = 12,        # 周期は12とする
        Q = NA
      ) +
      SSMarima(
        ar = c(0, 0),       # 2次のAR成分
        d = 0,
        Q = 0
      ) + distance + kuroshio_a + winter_a, # 外生変数
    data = ts_data
  )
  
  # optimに渡す前にパラメータをexpしたりartransformしたり、変換する
  # ほぼbuild_ssmと同じだが、パラメータだけ変更されている
  update_func <- function(pars, model) {
    model <- SSModel(
      H = exp(pars[6]),
      temp ~
        SSMtrend(degree = 2,
                 Q = c(list(0), list(exp(pars[1])))) +
        SSMseasonal(
          sea.type = "dummy",
          period = 12,
          Q = exp(pars[2])
        ) +
        SSMarima(
          ar = artransform(pars[3:4]),
          d = 0,
          Q = exp(pars[5])
        ) + distance + kuroshio_a + winter_a,
      data = ts_data
    )
  }
  
  # 最適化その1。まずはNelder-Mead法を用いて暫定的なパラメータを推定
  fit_ssm_bef <- fitSSM(
    build_ssm,
    inits = c(-17,-30, 0.5, 0, -1, -3), # パラメータの初期値(任意)
    update_func,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-16)
  )
  
  # 最適化その2。先ほどの結果を初期値に使ってもう一度最適化する
  fit_ssm <- fitSSM(
    build_ssm,
    inits = fit_ssm_bef$optim.out$par,
    update_func,
    method = "BFGS",
    control = list(maxit = 5000, reltol = 1e-16)
  )
  
  # フィルタリングとスムージング
  result_ssm <- KFS(
    fit_ssm$model,
    filtering = c("state", "mean"),
    smoothing = c("state", "mean", "disturbance")
  )
  
  # 結果の出力
  return(list(fit_ssm, result_ssm))
  
}

# 黒潮流路未使用
make_ssm_without <- function(ts_data) {
  # モデルの構造を決める
  build_ssm <- SSModel(
    H = NA,
    temp ~
      SSMtrend(degree = 2,                  # 平滑化トレンドモデル
               Q = c(list(0), list(NA))) +
      SSMseasonal(
        sea.type = "dummy", # ダミー変数を利用した季節成分
        period = 12,        # 周期は12とする
        Q = NA
      ) +
      SSMarima(
        ar = c(0, 0),       # 2次のAR成分
        d = 0,
        Q = 0
      ),
    data = ts_data
  )
  
  # optimに渡す前にパラメータをexpしたりartransformしたり、変換する
  # ほぼbuild_ssmと同じだが、パラメータだけ変更されている
  update_func <- function(pars, model) {
    model <- SSModel(
      H = exp(pars[6]),
      temp ~
        SSMtrend(degree = 2,                  # 平滑化トレンドモデル
                 Q = c(list(0), list(exp(pars[1])))) +
        SSMseasonal(
          sea.type = "dummy",
          period = 12,
          Q = exp(pars[2])
        ) +
        SSMarima(
          ar = artransform(pars[3:4]),
          d = 0,
          Q = exp(pars[5])
        ),
      data = ts_data
    )
  }
  
  # 最適化その1。まずはNelder-Mead法を用いて暫定的なパラメータを推定
  fit_ssm_bef <- fitSSM(
    build_ssm,
    inits = c(-17,-30, 0.5, 0, -1, -3), # パラメータの初期値(任意)
    update_func,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-16)
  )
  
  # 最適化その2。先ほどの結果を初期値に使ってもう一度最適化する
  fit_ssm <- fitSSM(
    build_ssm,
    inits = fit_ssm_bef$optim.out$par,
    update_func,
    method = "BFGS",
    control = list(maxit = 5000, reltol = 1e-16)
  )
  
  # フィルタリングとスムージング
  result_ssm <- KFS(
    fit_ssm$model,
    filtering = c("state", "mean"),
    smoothing = c("state", "mean", "disturbance")
  )
  
  # 結果の出力
  return(list(fit_ssm, result_ssm))
  
}


# モデル化 --------------------------------------------------------------------

# 黒潮流路のデータ利用
list_kuroshio <- make_ssm_kuroshio(train)

fit_kuroshio    <- list_kuroshio[[1]]
result_kuroshio <- list_kuroshio[[2]]


# 黒潮流路未使用
list_without <- make_ssm_without(train)

fit_without    <- list_without[[1]]
result_without <- list_without[[2]]


# 推定結果 --------------------------------------------------------------------

# 平滑化推定量
result_kuroshio$alphahat

# 係数
result_kuroshio$alphahat[1, 1:3]

# 係数の95%信頼区間
res <- confint(result_kuroshio, level = 0.95)
res[c("distance", "kuroshio_a", "winter_a")] %>% lapply(head, n = 1)

# 状態の可視化
autoplot(
  result_kuroshio$alphahat[, c("level", "slope", "sea_dummy1", "arima1")],
  facets = TRUE
)


# 予測 ----------------------------------------------------------------------

# 簡単な方法
predict(result_without$model, n.ahead = 14)

# newdataの利用

# 推定されたパラメータ
pars <- fit_without$optim.out$par

# newdataを用いた予測
pred_without <- predict(
  result_without$model, 
  newdata = SSModel(
    H = exp(pars[6]),
    rep(NA, nrow(validate)) ~
      SSMtrend(degree = 2,
               Q = c(list(0), list(exp(pars[1])))) +
      SSMseasonal(sea.type = "dummy",
                  period = 12,
                  Q = exp(pars[2])) +
      SSMarima(ar = artransform(pars[3:4]),
               d = 0,
               Q = exp(pars[5])),
    data = validate
  )
)

# 結果はまったく等しい
all(pred_without == predict(result_without$model, n.ahead = 14))

# なお、validateの代わりに、以下のtargetを用いても、予測結果は変わらない
# 精度の検証ではなく、本来的な将来予測を行う場合は、水温がNAとなっている方が自然であろう
# target <- validate
# target[, "temp"] <- NA




# 外生変数があるモデルの予測
# n.aheadは使えない
# predict(result_kuroshio$model, n.ahead = 14)と実行するとエラーになる

# newdataの利用
pars <- fit_kuroshio$optim.out$par

pred_kuroshio <- predict(
  result_kuroshio$model, 
  newdata = SSModel(
    H = exp(pars[6]),
    rep(NA, nrow(validate)) ~
      SSMtrend(degree = 2,
               Q = c(list(0), list(exp(pars[1])))) +
      SSMseasonal(sea.type = "dummy",
                  period = 12,
                  Q = exp(pars[2])) +
      SSMarima(ar = artransform(pars[3:4]),
               d = 0,
               Q = exp(pars[5])
      ) + distance + kuroshio_a + winter_a,
    data = validate
  )
)

# 参考
?predict.SSModel

# 予測値
pred_kuroshio
pred_without

# 予測精度
accuracy(pred_kuroshio, validate[, "temp"])
accuracy(pred_without,  validate[, "temp"])


