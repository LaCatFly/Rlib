################# utility session ##############

cap <- function (x, thres) {
  cap.tmp = na.omit(x[x > thres])
  dn = len(cap.tmp)
  if (dn != 0) {x[x > thres & !is.na(x)] = rep(thres, dn)}
  x 
}

consec.chg <- function(data, positive = T, thres = 0) {
  data.diff = data/mlag(data) - 1
  if(positive) dir = data.diff > thres else dir = data.diff < thres
  temp = cumsum(iif(dir, 1, 0))
  temp = temp - ifna.prev(iif(dir, NA, coredata(temp)))
  temp[temp == 0] = NA
  temp
}

normalize.median.mad <- function (x) { (x - median(x, na.rm = T))/mad(x, na.rm = T)}

expo.w <- function (nr, lambda = 0.994, order = TRUE){
  ## gen expotentially smoothing factor
  scal <- 1 - lambda
  tpR1.2 <- 0
  for (i in 1:nr){tpR1.2[i] <- lambda^i }
  
  expo.w <- tpR1.2/sum(tpR1.2)
  expo.w <- sort(expo.w, decreasing = order)
  expo.w
  
}

flr <- function (x, thres) {
  flr.tmp = na.omit(x[x < thres])
  dn = len(flr.tmp)
  if (dn != 0) {x[x < thres & !is.na(x)] = rep(thres, dn)}
  x
}

findOL <- function(xrow) {
  tmp = na.omit(xrow)
  mdn <- median(tmp)
  devs <- abs(tmp - mdn)
  return(which.max(devs))
}

find.sign <- function (x) { 
  tmp = x
  tmp[tmp > 0] = 1
  tmp[tmp == 0] = 0
  tmp[tmp < 0] = - 1
  tmp
}


ifna.zero <- function (data) {
  tmp  = is.nan(data) | is.na(data) 
  data[tmp] = 0 
  data
}

is.even <- function(x) {x %% 2 == 0}


re.order <- function (order.template, order.target) {
  tmp = order.target[match(order.template, order.target) ]
  tmp
}


rel.tan <- function(x) {
  tmp.err = is.null(x) | is.nan(x) | is.na(x) | is.infinite(x) | x == 0
  tmp = x[!tmp.err]
  if (len(tmp) == 0) {b = NA}  else { coeff = line(1:len(tmp), tmp)$coefficients
  b  = coeff[2] ; a  = coeff[1]  
  }
  return(b/a)
}

rel.tan2 <- function(x) {
  tmp.err = is.null(x) | is.nan(x) | is.na(x) | is.infinite(x) | x == 0
  tmp = x[!tmp.err]
  if (len(tmp) == 0) {b = NA} else { coeff = line(1:len(tmp), tmp)$coefficients; b  = coeff[2]}
  b
}

rel.tan.ols <- function(x) {
  tmp.err = is.null(x) | is.nan(x) | is.na(x) | is.infinite(x) | x == 0
  tmp = x[!tmp.err]
  if (len(tmp) == 0) {b = NA} else { coeff = ols((1:len(tmp)), tmp)$coefficients; b  = coeff[1]}
  b
}

rel.tan.func <- function(x, func, out.num) {
  tmp.err = is.null(x) | is.nan(x) | is.na(x) | is.infinite(x) | x == 0
  tmp = x[!tmp.err]
  tmp.func = match.fun(func)
  if (len(tmp) == 0) {b = NA} else { coeff = tmp.func(tmp , 1:len(tmp))$coefficients;}
  coeff[out.num]
}


ts.func.apply <- function (ts, func, col.row = 2, ...) 
{
  stopifnot(is.xts(ts))
  tmp.func = match.fun(func)
  tmp<-apply(X=ts, MARGIN = col.row, FUN=tmp.func, ...)
  if(!is.xts(tmp)) {tmp <- xts(tmp, order.by=index(ts))}
  ts.func.apply <- tmp
  return(ts.func.apply)
}

z.scores <- function(x) {
  (x- mean(x, na.rm=T)) / sd(x, na.rm=T)
}

z.scores.md <- function(x) {
  (x- median(x, na.rm=T)) / sd(x, na.rm=T)
}

rebase <- function(x){
  bs = first(x)
  rebase= apply(x, MARGIN =1, function(x) x/bs)
  return(t(rebase))
}


# ---------------------------- Portfolio Related ---------------------------------
  
fa.tptbm <- function (fac, gp.num = 3) {  
  tpf = seq(from = 1, to = len(ls(fac)), by = gp.num)
  bpf = seq(from = gp.num, to = len(ls(fac)), by = gp.num)
  tmp.fact = matrix(rep(fac[[1]]$ret, len(tpf)), ncol = len(tpf)) *NA
  colnames(tmp.fact) = paste0('P-',tpf,'-',bpf)
    
  for (i in 1:len(tpf )) {tmp.fact[, i] = fac[[tpf[i]]]$ret - fac[[bpf[i]]]$ret }
  tmp.fact = make.xts(rowSums(tmp.fact, na.rm = T)/len(tpf), order.by = index(fac[[1]]$ret))
  tmp.fact
} 


fa.tptbm <- function (fac, gp.num = 3) {  
  tpf = seq(from = 1, to = len(ls(fac)), by = gp.num)
  bpf = seq(from = gp.num, to = len(ls(fac)), by = gp.num)
  tmp.fact = matrix(rep(fac[[1]]$ret, len(tpf)), ncol = len(tpf)) *NA
  colnames(tmp.fact) = paste0('P-',tpf,'-',bpf)
  
  for (i in 1:len(tpf )) {tmp.fact[, i] = fac[[tpf[i]]]$ret - fac[[bpf[i]]]$ret }
  tmp.fact = make.xts(rowSums(tmp.fact, na.rm = T)/len(tpf), order.by = index(fac[[1]]$ret))
  tmp.fact
}

port.wgt.helper <- function(hist.ret, weight.algo.func, 
                            exp.ret = F, annual.factor = 252, 
                            l.b = 0, u.b = 0.1, 
                            cov.method = 'Ledoit' ) {
  library('tawny')
  if (exp.ret != F) { ia$expected.return = exp.ret}
  ia = create.historical.ia(hist.ret, annual.factor)
  s0 = apply(coredata(hist.ret), 2, sd)
  
  constraints = new.constraints(ia$n, lb = l.b, ub = u.b)
  constraints = add.constraints(rep(1, ia$n), 1, type = '=', constraints)        
  
  if (cov.method == 'kendall') {ia$cov = tawny::cov_sample(coredata(hist.ret))} 
  else {ia$cov = tawny::cov_shrink(coredata(hist.ret))}
  fc = match.fun(weight.algo.func)
  tmp = fc(ia, constraints)
  port.wgt.helper =  as.matrix(round(tmp , digits = 4))
  rownames(port.wgt.helper) = colnames(hist.ret)
  colnames(port.wgt.helper) = 'wght'
  port.wgt.helper
}

price.env.update <- function (env, sdate, to.Date, str = 'open,high,low,close,adjfactor') {
  bt.env.clean(env)
  tks = ls(env)
  stmp = paste0('::', as.Date(sdate) - 1)
  for (i in tks) {
    tmp = w.xts(w.wsd(i, str, sdate, to.Date))
    data.wd.price.only[[i]] = rbind.xts(env[[i]][stmp], tmp)
    print(i)}
  save(data.wd.price.only, file = paste0(getwd(), paste0('/data-price-only',as.Date(to.Date),'.Rdata')))
}

attr.env.update <- function (env, sdate, to.Date, update.str){
  bt.env.clean(env)
  tks = ls(env)
  if (ncol(env[[tks[1]]]) !=  len(spl(update.str)) ) {stop('update.str is not equal length')}
  str = paste0('::', sdate)
  stmp = paste0('::', as.Date(sdate) - 1)

  for (i in tks) {
    tmp = w.xts(w.wsd(i, update.str, sdate, to.Date, "unit=1;traderType=1;shareType=1"))
    data.wd.attr[[i]] = rbind.xts(env[[i]][stmp], tmp)
    print(i)
  }
  save(data.wd.attr, file = paste0(getwd(), paste0('/data-attr',as.Date(to.Date),'.Rdata')))
}

env.price.adj <- function (x.env) {
  tks = ls(x.env)
  p.adj = new.env()
  for (i in tks) {
    tmp = x.env[[i]]
    cur.adj.factor = (rep(coredata(tmp[, 5]), 4) %>% matrix( ncol = 4))/(last(tmp)[ ,'ADJFACTOR'] %>% as.numeric) 
    date.idx = index(tmp)
    p.adj[[i]] = make.xts(coredata(tmp[ , -5]) * cur.adj.factor, order.by = date.idx)
  }
  return(p.adj)
}

env.cbind <- function (clean.env1, clean.env2) {
  ls1 = ls(clean.env1)
  ls2 = ls(clean.env2)
  if (!identical(ls1, ls2)) {stop ('env index are not at the same length')}
  date.idx1 = clean.env1[[ls1[1]]] %>% index %>% as.Date
  date.idx2 = clean.env2[[ls2[1]]] %>% index %>% as.Date
  if (any(date.idx1 != date.idx2)) {stop ('date index are not identical')}
  for (i in ls1) {clean.env1[[i]] = cbind.xts(clean.env1[[i]], clean.env2[[i]])}
  return(clean.env1)
}


env.merge <- function (clean.env1, clean.env2) {
  ls1 = ls(clean.env1)
  ls2 = ls(clean.env2)
  if (!identical(ls1, ls2)) {stop ('env index are not at the same length')}
  
  for (i in ls1) {clean.env1[[i]] = cbind.xts(clean.env1[[i]], clean.env2[[i]])}
  return(clean.env1)
}

env.check <- function (env.p  =  data.wd.price.only, env.attr = data.wd.attr, update.str = est) {
  tks.p = ls(env.p)
  tks.attr = ls(env.attr)
  print(paste0('p.env is the same as attr.env? ', identical(tks.p, tks.attr )) )
  
  tmp.p = c(1: len(tks.p) ) * NA 
  for (j in 1:len(tks.p) ) { tmp.p[j] = ncol(env.p[[tks.p[j]]]) }
  print(paste0('all tks in p.env are the same col? ', any(tmp.p == tmp.p[1]), ' @ ', tmp.p[1]) )
  
  tmp.attr = c(1: len(tks.attr) ) * NA 
  for (j in 1:len(tks.attr) ) { tmp.attr[j] = ncol(env.attr[[tks.attr[j]]]) }
  print(paste0('all tks in attr.env are the same col? ', any(tmp.attr == tmp.attr[1]), ' @ ', tmp.attr[1]) )
  print(paste0('all tks in attr.env are the same dim as update.str? ', any(tmp.attr == len(spl(update.str)))) ) 
  print(paste0('tks colnames are identical to update.str? ', identical(colnames(env.attr[[tks.attr[j]]]), toupper(spl(update.str)))) )
}

env.ifna.prev <- function (env, to.replace) {
  ls1 = ls(env)
  for (i in ls1) {env[[i]] = na.locf(env[[i]])}
  return(env)
}

env.rmColto <- function (env, ncol.name) {
  ls1 = ls(env)
  toEnv = new.env()
  col.str = colnames(env[[ls1[1]]])
  colrm = match(toupper(ncol.name), col.str)
  for (i in ls1) {
    toEnv[[i]] = env[[i]][, ncol.name]
    env[[i]] = env[[i]][, -colrm]
  }  
  return(toEnv)
}


bt.chart <- function (x, benchmark, time.windows) {
  tmp.bind = NULL
  for (i in ls(x)) {tmp.bind = merge.xts(tmp.bind, x[[i]]$ret)}
  tmp.bind = merge.xts(tmp.bind, benchmark)[time.windows]
  tmp.bind
}


bt.ts.grp.split.l3 <- function (factor.list, reval.date, factor.order = c(1,2,3), ex.list = ivol.tmp, n.grp = 2) {
  f.n = ls(factor.list)
  if (len(f.n) != len(factor.order)) { stop ('order list is not equal to factor list')}
  sd.factor = list(); tmp.grp = list(); tmp.grp.l2 = list(); tmp.grp.l3 = list();
  ex.tmp = which(ex.list[reval.date, ] == 0 | is.na(ex.list[reval.date, ]) | is.nan(ex.list[reval.date, ]) )
  for (k in f.n ) {
    sd.factor[[k]] = t(factor.list[[k]][reval.date, colnames(factor.list[[k]])[-ex.tmp]])
  }
  
  tmp.grp = bt.grp.split(rownames(sd.factor[[factor.order[1]]]), n.grp = n.grp, grp.by = sd.factor[[factor.order[1]]], label = 'grp')
  for (j in ls(tmp.grp)) {
    f1.byref = sd.factor[[factor.order[2]]][match(tmp.grp[[j]], rownames(sd.factor[[factor.order[2]]])), ] %>% transform
    tmp.grp.l2[[j]] = bt.grp.split(tmp.grp[[j]], n.grp = n.grp, grp.by = f1.byref, label = j)
    for (g in ls(tmp.grp.l2[[j]])) {
      f2.byref = sd.factor[[factor.order[3]]][match(tmp.grp.l2[[j]][[g]], rownames(sd.factor[[factor.order[3]]])), ] %>% transform
      tmp.grp.l3[[j]][[g]] = bt.grp.split(tmp.grp.l2[[j]][[g]], n.grp = n.grp, grp.by = f2.byref, label = g)
    }
  }
  return( unlist(tmp.grp.l3, recursive = F, use.names = F) %>% unlist(recursive = F))
}

# ---------------------------- Trading related ----------------------------------------

gen.trd.score <- function (pob.hi, pob.lo, pob.cl, pob.volu){
  ohc.ob = cbind(pob.hi, pob.cl, pob.volu)
  nrow = dim(ohc.ob)[1]
  attr.volu = RollingZscore( EMA(pob.volu, 3), window = 100, na_rm = T) %>% pnorm
  volu.p.cnsec = consec.chg(attr.volu, positive = F)
  volu.n.cnsec = consec.chg(attr.volu)
  
  vu.m = RollingMean(attr.volu, window = 30, na_rm = T)
  vu.std = RollingStd(attr.volu, window = 150, na_rm = T)
  
  vu.p.ceil = vu.m + 2*vu.std
  vu.p.ex = vu.m + 0.2*vu.std
  vu.n.ex = vu.m - 0.3*vu.std
  vu.n.flr = vu.m - 0.7*vu.std
  
  volu.l = iif(attr.volu > vu.p.ex & attr.volu < vu.p.ceil, 1, NA) %>% sig.clean(0, 2)
  volu.s = iif(attr.volu < vu.n.ex & attr.volu > vu.n.flr, -1, NA) %>% sig.clean(0, 3)
  
  mon.lt = ifelse(EMA(pob.cl, 5)> EMA(pob.cl, 10) & EMA(pob.cl, 10)> EMA(pob.cl, 20), 1, NA) %>% sig.clean(0, 1)
  mon.st = (ifelse(EMA(pob.cl, 3) > EMA(pob.cl,5) , 1, NA)) %>% sig.clean(0, 1)
  mon.rr = ifelse(EMA(pob.cl, 10) < EMA(pob.cl, 20), -1, NA) %>% sig.clean( 0, 1)
  
  vu.trnd.p = iif(volu.p.cnsec > 1 & (EMA(pob.cl, 3) > EMA(pob.cl,5)) , 1, NA) %>% sig.clean(0, 1)
  vu.trnd.n = iif(volu.n.cnsec > 2 & (EMA(pob.cl, 5) < EMA(pob.cl,20)) , -1, NA) %>% sig.clean(0, 1)
  
  mfi5 = MFI(ohc.ob, volume = pob.volu, n = 8)
  mfi3 = MFI(ohc.ob, volume = pob.volu, n = 3)
  mfi3.l2 = (mfi3/mlag(mfi3, 2) - 1)
  mfi.l = iif(mfi3.l2 < 0 & mfi5  > 75, 1, NA) %>% sig.clean(0, 2)
  
  sig.l = cbind(vu.trnd.p, volu.l, mon.lt , mon.st, mfi.l)
  sig.s = cbind(vu.trnd.n, volu.s, mon.rr)
  sig.all = cbind(sig.l, sig.s)
  wgh.tmp = c(0.2, 0.15, 0.15, 0.15, 0.15, 0.3, 0.2, 0.2)
  wght = rep(wgh.tmp/sum(wgh.tmp), nrow) %>% matrix(ncol = ncol(sig.all), byrow = TRUE)
  
  scores = (sig.all * wght) %>% apply(MARGIN = 1, sum)
  return(scores)
}



################# Wind related ##############
w.xts <- function(data, na.prev = F) {
  if (data$ErrorCode!= 0) {stop ("data$ErrorCode")}
  if(is.null(data$ErrorCode) | is.null(data$Data)) { stop('not wind object') }
  core.data = data$Data[, -1]
  if (na.prev) { 
      if (is.null(dim(core.data))) { core.data = ifna.prev (core.data)} else { 
          core.data = apply(core.data, MARGIN = 2, ifna.prev)  }
    }  
  tmp = make.xts(core.data , order.by = data$Data[, 1] )
  if (!is.null(data$Code) & dim(tmp)[2] == 1) { colnames(tmp) = data$Code}
  tmp
}


w.xts.convert <- function(data) {
  tmp=as.Date(data[, 1], "%y/%m/%d")
  ts.data = make.xts(data[ ,-1], order.by=tmp)
  return(ts.data) 
}

w.format.convert <- function(x) {
  if(dim(x$Data)[2] == 2 ) { 
      tmp.name = x$Code
      colnames(x$Data)[2] = tmp.name
      x$Code = 'CLOSE' 
      }
  return(x)
}

w.as.data <- function(wdata , origin = "1900-01-01") {
  if(wdata$ErrorCode != 0) {stop ('wind data err')}
  for (i in 2:(dim(wdata$Data)[2])) {
    tmp= as.Date(wdata$Data[ , i], origin = origin ) 
    wdata$Data[, i] = tmp
  }
  wdata$Data
}

w.datetoStr <- function(x) {
  tmp = x
  if (class(x) != 'Date') {tmp = as.Date(x)}
  tmp = gsub("-", replacement = "", as.character(tmp))
  tmp
} 

w.f.rpt.date <- function(x) {
  tmp = x
  if (class(x) != 'Date') {tmp = as.Date(x)}
  date.dff =  tmp - rpt.date.yr.db
  w.f.rpt.date = rpt.date.yr.db[min(which(date.dff > 0))]
  return(w.f.rpt.date)
}

w.sort.rpt.date <- function (rpt.data, def.avg = 30) {
  tmp.index = index(rpt.data)
  tmp = apply( rpt.data , MARGIN = 2, function(x) as.Date(x) - tmp.index) 
  for (i in 1:dim(tmp)[2]) { tmp[ ,i][tmp[ , i] < 0] = NA }
  avg = apply(tmp, MARGIN = 2, function(x) mean(x, na.rm = TRUE) ) %>% round(digits = 0)
  for (i in 1:dim(tmp)[2]) { 
    na.ct = len(which(!is.na(tmp[ , i])))
    if (na.ct < 2 ) {tmp[, i] = def.avg } else {tmp[is.na(tmp[ , i]), i] = avg[i] }
  }
  w.sort.rpt.date = make.xts(tmp, order.by = tmp.index)
  w.sort.rpt.date
}

w.qfix.date <- function(date.list, ref.data, interval = 50) {
  tmp = (ref.data - date.list - interval) 
  w.qfix.date = d.list[which.min(tmp[tmp >0])]
  w.qfix.date
}

w.data.list <- function(x) {
  tmp = list()
  for (i in colnames(x)) {tmp[[as.character(i)]] = x[,i ] }
  return(tmp)
}


w.append.ind <- function(ticker, trdeDate.yyyy.mm.dd) {
  
  d.tmp = paste0('tradeDate=', trdeDate.yyyy.mm.dd)
  ind.map = w.wss(unique(ticker),'industry_sw,indexcode_sw','industryType=1',d.tmp)$Data
  ind.tmp = ind.map[match(ticker, ind.map$CODE), ]
  ind.tmp
}

wd.prep <- function(env) {
  sym = ls(env)
  ind = index(env[[sym[1]]])
  tk =  colnames(env[[sym[[1]]]])
  wd = new.env()
  
  for (i in tk) {
    tmp = list()
    for (icol in sym) { tmp[[icol]] = env[[icol]][ , i] }  
    tk.tmp = make.xts(abind(tmp, along = 2), order.by = ind)
    colnames(tk.tmp) = toupper(sym)
    wd[[i]] = tk.tmp
  }
  wd
}

########## BBG related ########## 
bbg.data.pre <-function (bbg.data, env, col.name = 'Close') {
  if (class(bbg.data)!= 'list') {stop ('input data is not list')}
  tmp = lapply(bbg.data, function (x) if (!is.null(x)) {as.xts(x[, -1], order.by = x[, 1] )} )
  for(i in ls(tmp)) {
    colnames(tmp[[i]]) = col.name
    env[[i]] = tmp[[i]]
  } 
}



########### Data clean ########### 

dc.iflarge.prev <- function(data, threshold =5) {
  tmp = z.scores(data)
  sm = abs(tmp) < threshold
  out = which(abs(tmp) > threshold)
  if (len(out) >0) {cat ('outliners are identified \n', 
         colnames(data),":",tmp[out], '\n',
         'position :', out, '\n')}
  sm[1] = T
  st = first(which(!is.na(data)))
  data[cummax( (st:len(data)) * sm[st:len(sm)])]
}



dc.shift <- function(data, min.ratio = 2.5) {
  #used for SS M&A case
  ticker.data = data
  nperiods = nrow(ticker.data)
  price = ticker.data
  ratio = as.vector((price)/mlag(price))
  index = which(ratio > min.ratio)
  if(len(index) > 0) {
    for(j1 in index) {
      cat('Abnormal price found for', colnames(c.data[, i]), 
          format(index(ticker.data)[j1],'%d-%b-%Y'),'Ratio :', 
          round(ratio[j1],1),"position: ", j1,":", i, '\n')
      ticker.data[j1:nperiods] = ticker.data[j1:nperiods] / ratio[j1]
    }
  }  
  ratio = as.vector(mlag(price)/(price))
  index = which(ratio > min.ratio)
  if(len(index) > 0) {
    for(j2 in index) {
      cat('Abnormal price found for',  colnames(c.data[, i]), 
          format(index(ticker.data)[j2],'%d-%b-%Y'),'Inverse Ratio :', 
          round(ratio[j2],1),"position: ", j2,":", i, '\n')
      ticker.data[j2:nperiods] = ticker.data[j2:nperiods] * ratio[j2]
    }
  }
  return(ticker.data)
}


dc.w.xts <- function(data, cut.off = 5){ 
  c.data = data[rowSums(is.na(data)) == 0, ]
  for (i in 1: ncol(c.data)) {
    tmp = c.data[, i]
    tmp = dc.iflarge.prev(tmp, cut.off)
    c.data[, i] = tmp
  }
  dc.w.xts  = c.data
}

dc.full.normalize <- function (x, qt = c(0.97, 0.03)) {
  tmp = ifna(x, NA)
  tmp = (tmp - median(x, na.rm = T))/mad(x, na.rm =T)
  q = quantile(tmp, prob = qt, na.rm = T)
  tmp[which(tmp > q[1])] = q[1] ;  tmp[which(tmp < q[2])] = q[2]
  tmp
}

dc.abs.sqrt <- function(x, tar = 2.3, diff.thres = 0.5) {
  xdata = as.matrix(x)
  infc = is.infinite(xdata)
  if (any(infc)) { xdata[which(infc)] = NA }
  if (sum(is.nan(xdata), is.na(xdata)) == len(xdata)) {tmp = x } 
  else {
          tmp = abs(xdata); s = xdata/tmp
          extr = max(xdata, na.rm = T) - abs(min(xdata, na.rm = T))
          tail.tmp = max(max(xdata, na.rm = T), abs(min(xdata, na.rm = T)) )
          if (abs(extr) < diff.thres | tail.tmp > tar) {tmp = (tmp)^dc.pow(tail.tmp, tar) * s} else {tmp  = xdata} }
  tmp
}

dc.pow <- function (x, target = 2.75) {
  if (!is.infinite(x) ) {tmp = solve(log10(x), log10(target))} else { tmp = NA}
  tmp
}

dc.rm.static <- function (x) {
  x[ -unique(match(x, x))] = NA
  x
}


dc.bump.mid <- function(x, pow = 0.6) {
  tmp = 1/(abs(x)^pow) 
  tmp
}
################# Calc method sessions###################
ri.csfb<- function(data, 
                    w.ret = 15,  #short term rolling return windows
                    w.vol = 25, 
                    l.vol.win = 250, 
                    weights = rep(1, dim(data)[2]),
                    log.flag = FALSE,  rev.flag = FALSE) 
{

  h.ret10 = bt.apply.matrix(data, ROC, n = w.ret, type =  'discrete') %>% na.locf
  h.ret1 = bt.apply.matrix(data, ROC, type =  'discrete') %>% na.locf
  h.vol = bt.apply.matrix(h.ret1, runSD, n = w.vol)
  l.vol = bt.apply.matrix(h.ret1, runSD, n= l.vol.win )
  
  z.ret = (h.ret10 - l.vol)/(l.vol)
  p.vol = h.vol/l.vol 
  if(log.flag) {p.vol = log(h.vol/l.vol)} 
 
  
  b1 = h.ret1[, 1:3] * NA
  colnames(b1) = c('Value', 'Std.Err', 't-value')
  
  for (i in 1:dim(h.ret1)[1]) {
    vol = t(p.vol[i, ])
    rt = z.ret[i, ] * weights
    
    vol1 = t(p.vol[i, which(!is.na(vol))])
    rt1 = t(rt[, which(!is.na(vol))])
    
    if (rev.flag) {tmp = vol1; vol1 = rt1; rt1 = tmp}
    
    msg = try(rlm(vol1, rt1,  psi = psi.bisquare,  wt.method = 'inv.var'), silent = T)
    if (class(msg)[1] == 'try-error') {
      warning(i, msg, '\n')
    } else {b1[i, ] = summary(msg)$coefficients}
  }

  ri = list(); ri$rt = z.ret; ri$vol = p.vol; ri$b1 = b1
  return(ri)
}



ri.fti <- function(data, look.back = 90) {
  nperiods = nrow(data)
  turbulence = data[,1] * NA
  ret = coredata(data / mlag(data) - 1)
  
  for( i in (look.back+1) : nperiods ) {
    temp = ret[(i - look.back + 1):(i-1), ]
    
    turbulence[i] = mahalanobis(ret[i,], colMeans(temp), cov(temp))
    if( i %% 200 == 0) cat(i, 'out of', nperiods, '\n')
  } 
  ri.fti = turbulence
}



bt.matching.find2 <- function
(
  data,
  n.query=90,
  n.reference=252*10,
  n.match=10,
  normalize.fn = normalize.mean.sd,
  dist.fn = dist.euclidean,
  plot=TRUE,
  plot.dist=FALSE,
  layout = NULL,
  main = NULL
)
{
  data = last(data, n.reference)
  reference = coredata(data)
  n = len(reference)
  query = reference[(n - n.query + 1):n]
  reference = reference[1:(n - n.query)]
  main = paste(main, join(format(range(index(data)[(n - n.query + 1):n]), '%d%b%Y'), ' - '))
  n.query = len(query)
  n.reference = len(reference)
  dist.fn.name = ''
  if(is.character(dist.fn)) {
    dist.fn.name = paste('with',dist.fn)
    dist.fn = get(dist.fn)
  }
  dist = rep(NA, n.reference)
  query.normalized = match.fun(normalize.fn)(query)
  for( i in n.query : n.reference ) {
    window = reference[ (i - n.query + 1) : i]
    window.normalized = match.fun(normalize.fn)(window)
    dist[i] = match.fun(dist.fn)(rbind(query.normalized, window.normalized))
    if( i %% 100 == 0) cat(i, '\n')
  }
  min.index = c()
  temp = dist
  temp[ temp > mean(dist, na.rm=T) ] = NA
  for(i in 1:n.match) {
    if(any(!is.na(temp))) {
      index = which.min(temp)
      min.index[i] = index
      temp[max(0,index - 2*n.query) : min(n.reference,(index + n.query))] = NA
    }
  }
  n.match = len(min.index)
  if(plot) {
    dates = index(data)[1:len(dist)]
    if(is.null(layout)) {
      if(plot.dist) layout(1:2) else layout(1)
    }
    par(mar=c(2, 4, 2, 2))
    if(plot.dist) {
      plot(dates, dist, type='l',col='gray', main=paste('Top Historical Matches for', main, dist.fn.name), ylab='Distance', xlab='')
      abline(h = mean(dist, na.rm=T), col='darkgray', lwd=2)
      points(dates[min.index], dist[min.index], pch=22, col='red', bg='red')
      text(dates[min.index], dist[min.index], 1:n.match, adj=c(1,1), col='black',xpd=TRUE)
    }
    plota(data, type='l', col='gray', LeftMargin = 1,
          main=iif(!plot.dist, paste('Top Historical Matches for', main), NULL)
    )
    plota.lines(last(data,n.query), col='blue')
    for(i in 1:n.match) {
      plota.lines(data[(min.index[i]-n.query + 1):min.index[i]], col='red')
    }
    text(index4xts(data)[min.index - n.query/2], reference[min.index - n.query/2], 1:n.match,
         adj=c(1,-1), col='black',xpd=TRUE)
    plota.legend(paste('Pattern: ', main, ',Match Number'),'blue,red')
  }
  return(list(min.index=min.index, dist=dist[min.index], query=query, reference=reference, dates = index(data), main = main))
}


# FAA Keller & Putten
FAA.scores <- function(prices, bench, wR = 1, wV = 0.5, wC = 0.8,  n.ret = 3 * 22, n.vol = 2 * 22, n.cor = 3 * 22) {
  if (dim(prices)[1] != dim(bench)[1]) { stop ('inputs dimension are not the same')}
  week.ends = endpoints(prices, 'weeks')
  week.ends = week.ends[week.ends > 0] 
  ri.wk = prices[week.ends,]
  ri.wk.ret = ri.wk/mlag(ri.wk) - 1
  ret = prices/mlag(prices) - 1
  intrv.wk = ivol.tmp[week.ends, ] - mlag(ivol.tmp[week.ends, ])
  
  ri = prices / mlag(prices, n.ret)
  rmax = apply(ri, MARGIN = 1, max, na.rm = T) #%>% transform
  m = log(ri/rmax) 
  
  
  vi = bt.apply.matrix(ret, runSD, n = n.vol)
  vmin = apply(vi, MARGIN = 1, function (x) min(x[x>0.0005], na.rm = T) ) 
  v = ifna(-sqrt(-log(vmin/vi)), NA)
  v = apply(v, MARGIN = 2, cap, 1)
  v = make.xts(v, order.by = index(prices))
  
  dl = dim(prices)
  cor.mat =  bt.apply.matrix(ret, function (x) RollingCorr(x, bench, n.cor, na_rm = T)) 
  
  # for (i in ( n.cor + 1):dl[1]) {
  #   hist = ret [(i - n.cor + 1):i, ]
  #   hcor = cor.shrink(hist)
  #   cor.idx = index(ret[i, ]) %>n% as.Date
  #   cor.mat[cor.idx, ] = apply(hcor, MARGIN = 2, function(x) (sum(x) - 1)/(dl[2]- 1) )
  # }
  
  ci = ifna(cor.mat, NA)
  cmin = apply(ci, MARGIN = 1, function (x) min((x + 0.0005), na.rm = T))
  c1 = log((cmin + 1)/(ci + 1))
  
  li = wR * m + wV * v  + wC * c1
  li
}


faa.scores2 <- function(prices, swing, bench.price, thres = 50, wR = 1, wV = 0.5, wC = 0.8,  n.ret = 12, n.vol = 8, n.cor = 12) {
  if (dim(prices)[1] != dim(bench.price)[1]) { stop ('inputs dimension are not the same')}
  week.ends = endpoints(prices, 'weeks')
  week.ends = week.ends[week.ends > 0] 
  ri.wk = prices[week.ends,]
  ri.wk.ret = ri.wk/mlag(ri.wk) - 1
  ret = prices/mlag(prices) - 1
  intrv.wk = ivol.tmp[week.ends, ] - mlag(ivol.tmp[week.ends, ])
  bch.wk.ret = bench.price[week.ends, ]/mlag(bench.price[week.ends, ]) - 1
  
  ri = ri.wk/ mlag(ri.wk, n.ret)
  rmax = apply(ri, MARGIN = 1, max, na.rm = T) #%>% transform
  m = log(ri/rmax) 
  
  vi = ri.wk * NA; ci = ri.wk * NA
  dn = dim(ri.wk)[2]
  for (i in 1:dn) {
    tmp.wk.ret = ri.wk.ret[ intrv.wk[, i] != 0,  i]
    if (len(tmp.wk.ret) < thres ) {tmp.sd = NA; tmp.cor = NA} 
    else {
      tmp.sd = RollingStd(tmp.wk.ret, n.vol, na_rm = T)
      tmp.cor = RollingCorr(tmp.wk.ret, bch.wk.ret[intrv.wk[, i] != 0, ], n.cor, na_rm = T)    
    }
    tmp.sd = merge.xts(tmp.sd, ri.wk.ret[, i])[, 1]
    tmp.cor = merge.xts(tmp.cor, ri.wk.ret[, i])[, 1]
    vi[, i] = tmp.sd %>% ifna(NA)
    ci[, i] = tmp.cor %>% ifna(NA)
    print(i)
  }
  
  vmin = apply(vi, MARGIN = 1, function (x) min(x[x>0.0005], na.rm = T) ) 
  v = ifna(-sqrt(-log(vmin/vi)), NA)
  v = apply(v, MARGIN = 2, cap, 1)
  v = make.xts(v, order.by = index(ri.wk))
  
  cmin = apply(ci, MARGIN = 1, function (x) min((x + 0.0001), na.rm = T))
  c1 = log((cmin + 1)/(ci + 1))
  
  li = wR * m + wV * v  + wC * c1
  li
}


#################### BT related ################
bt.grp.split.get <- function (x, grp.by, n.grp = 10, get.grp) {
  out.tmp = NULL 
  na.grp.by = as.matrix(grp.by)
  x1 = as.matrix(x)
  if (dim(x1)[1] != dim(na.grp.by)[1]) {stop ('length is not the same') }
  grp = list()
 
  na.test = is.na(na.grp.by) | is.nan(na.grp.by)
  if (any(na.test)) {
    ref.tmp = na.grp.by[- which(na.test)] 
    tar.obj = x1[- which(na.test) ]
    }  else {
        ref.tmp = na.grp.by
        tar.obj = x1
        }
  qt = quantile(ref.tmp, probs = (1:n.grp)/n.grp, na.rm = T)
  for (i in 1: len(qt)) {
    if (len(qt[i-1]) == 0) {tmp = tar.obj [ref.tmp <= qt[i] ] }  
      else {tmp= tar.obj [ref.tmp > qt[i-1] & ref.tmp <=qt[i] ] }
    g.name = paste0('group-', i) 
    grp[[g.name]] = tmp
  }

  for (j in get.grp) {
    out.tmp = c(out.tmp, grp[[j]])
  }
  out.tmp
}


bt.grp.split <- function (x, grp.by, n.grp = 3, label = 'group') {
  out.tmp = NULL 
  grp = list()
  na.grp.by = ifna(as.matrix(grp.by), NA)
  x1 = ifna(as.matrix(x), NA)

  na.test = is.na(na.grp.by) | is.nan(na.grp.by)
  if (any(na.test)) {
    ref.tmp = na.grp.by[- which(na.test)] 
    tar.obj = x1[- which(na.test) ]
  }  else {
    ref.tmp = na.grp.by
    tar.obj = x1
  }
  
  if (dim(x1)[1] != dim(na.grp.by)[1]) {stop ('length is not the same') }
  
  qt = quantile(ref.tmp, probs = (1:n.grp)/n.grp, na.rm = T)
  for (i in 1: len(qt)) {
    if (len(qt[i-1]) == 0) {tmp = tar.obj [ref.tmp <= qt[i] ] }  
    else {tmp= tar.obj [ref.tmp > qt[i-1] & ref.tmp <=qt[i] ] }
    grp[[paste0(label, '-', i) ]] = tmp
  }
  grp
}


bt.ts.grp.split.l3 <- function (factor.list, reval.date, factor.order = c(1,2,3), ex.list = ivol.tmp, n.grp = 2) {
  f.n = ls(factor.list)
  if (len(f.n) != len(factor.order)) { stop ('order list is not equal to factor list')}
  sd.factor = list(); tmp.grp = list(); tmp.grp.l2 = list(); tmp.grp.l3 = list();
  ex.tmp = which(ex.list[reval.date, ] == 0 | is.na(ex.list[reval.date, ]) | is.nan(ex.list[reval.date, ]) )
  for (k in f.n ) {
    sd.factor[[k]] = t(factor.list[[k]][reval.date, colnames(factor.list[[k]])[-ex.tmp]])
  }
  
  tmp.grp = bt.grp.split(rownames(sd.factor[[factor.order[1]]]), n.grp = n.grp, grp.by = sd.factor[[factor.order[1]]], label = 'grp')
  for (j in ls(tmp.grp)) {
    f1.byref = sd.factor[[factor.order[2]]][match(tmp.grp[[j]], rownames(sd.factor[[factor.order[2]]])), ] %>% transform
    tmp.grp.l2[[j]] = bt.grp.split(tmp.grp[[j]], n.grp = n.grp, grp.by = f1.byref, label = j)
    for (g in ls(tmp.grp.l2[[j]])) {
      f2.byref = sd.factor[[factor.order[3]]][match(tmp.grp.l2[[j]][[g]], rownames(sd.factor[[factor.order[3]]])), ] %>% transform
      tmp.grp.l3[[j]][[g]] = bt.grp.split(tmp.grp.l2[[j]][[g]], n.grp = n.grp, grp.by = f2.byref, label = g)
    }
  }
  return( unlist(tmp.grp.l3, recursive = F, use.names = F) %>% unlist(recursive = F))
}



bt.grp.fact.test <- function (alpha.obj, test.env, func = func, 
                              lp = 5, date.ind = m.ends[-1:-30],
                              susp = ivol.tmp, out.name = 'alpha', 
                              trade.flag = F, capital = 100000000) {
  
  mod = list(); fa = list()
  for (j in ls(alpha.obj)) {fa[[j]] = as.xts(apply(alpha.obj[[j]], MARGIN = 2, dc.full.normalize))}
  for (k in 1:lp) {
    test.env$weight[] = NA * bt.apply.matrix(test.env$prices, function (x) ifna(x, NA))
    for (i in 1:len(date.ind)) {
      mi = date.ind[i]
      re.date =  index(prices[mi])
      ref.tmp  = lapply(fa, function (x) x[re.date]) %>% abind(along =1) %>% t()
      tk.nm = rownames(ref.tmp)
      ref.tmp = apply(ref.tmp, MARGIN = 2, dc.abs.sqrt) %>% apply( MARGIN = 1, func ) %>% transform    
      #apply(ref.tmp, MARGIN = 2, dc.abs.sqrt) %>%
      rownames(ref.tmp) = tk.nm
      err.tmp = coredata(susp[re.date, ])
      
      ex.tmp = which(err.tmp == 0 | is.na(err.tmp) | is.nan(err.tmp) )
      tk.ex.nm = rownames(ref.tmp)[-ex.tmp]
      ref.tmp = ref.tmp[tk.ex.nm, ] %>% transform
      rownames(ref.tmp) = tk.ex.nm
      grp = bt.grp.split(rownames(ref.tmp) , grp.by =ref.tmp , n.grp = lp)
      
      if (dim(ref.tmp)[1] == 0) { test.env$weight[re.date, wght] = NA} 
      else {
        if (i - 1 == 0) { wuc = 0} else {wuc = na.locf(test.env$weight[paste0('::', re.date), last.ex.tmp])[re.date, ] %>% sum(na.rm = T)}
        wght = match(grp[[k]], colnames(test.env$weight))
        dn = len(wght)
        test.env$weight[re.date, -ex.tmp] = 0 
        test.env$weight[re.date, wght] = cap(rep((1-wuc)/dn, dn), 0.05)
        last.ex.tmp = ex.tmp 
      }
      pw = round(sum(test.env$weight[re.date, ], na.rm = T), 4)
      print(paste0(k, '-mod-', re.date, ';  tw = ', pw))
    }
    m.name = paste0(out.name, k)
    mod[[m.name]] = bt.run.share(test.env, clean.signal=F, capital = capital, commission =  0.002, trade.summary = trade.flag)
    colnames(mod[[m.name]]$ret) = m.name
  }
  return(mod)
}

bt.port.opt.helper <- function(hist.ret, weight.algo.func, 
                               exp.ret = F, annual.factor = 252, 
                               l.b = 0, u.b = 0.1, 
                               cov.method = 'Ledoit' ) {
  ia = create.historical.ia(hist.ret, annual.factor)
  if (exp.ret != F) { ia$expected.return = bt.port.exp.ret(ia, ret.adj = F)}
  s0 = apply(coredata(hist.ret), 2, sd)
  
  if (cov.method == 'kendall') 
  {
    ia$cov = cor(coredata(hist), use='complete.obs', method=cov.method) * (s0 %*% t(s0))
  } else {ia$cov = tawny::cov_shrink(coredata(hist.ret))}
  
  constraints = new.constraints(ia$n, lb = l.b, ub = u.b)
  constraints = add.constraints(rep(1, ia$n), 1, type = '=', constraints)        
  
  fnc = match.fun(weight.algo.func)
  tmp = fnc(ia, constraints)
  bt.port.opt.helper =  as.matrix(round(tmp , digits = 4))
  rownames(port.opt.helper) = colnames(hist.ret)
  colnames(port.opt.helper) = 'wght'
  bt.port.opt.helper
}

bt.port.exp.ret <- function(x, ret.adj = TRUE, rec.wind = 4) {
  tmp = sqrt(abs(1+x$expected.return)) - 1
  if (ret.adj == TRUE) { tmp.n = round(x$n/2, digits = 0); 
  tmp.rec = tail(x$hist.returns, rec.wind) %>% apply(MARGIN = 2, function(x) cumprod(x +1) - 1) %>% tail(1);
  tmp = tmp/x$n * (1 - tmp.n) + tmp.rec/x$n * tmp.n
  }
  tmp
}

bt.min.cov.calc.sel <- function (wk.ret.input, look.back.start.date, look.back.win = 150) {
  s.date = as.Date(look.back.start.date)
  date.len = paste0( s.date - look.back.win,'::', s.date)
  tmp.prices = wk.ret.input[date.len, ]
  dn = dim(tmp.prices)[1]
  nc =  apply(tmp.prices, MARGIN = 2, function (x) sum((x!= 0) == T ))
  tmp.prices.rm = tmp.prices[ ,which(nc == dn)]
  tmp.prices.rm
}

bt.ex.col <- function(env, col.name, sym.flag = T) {
  if (sym.flag == T) {sym = env$symbolnames} else { sym = ls(env) }
  if (is.null(sym) & sym.flag == T) { bt.prep(env); sym = env$symbolnames }
  tmp = list()
  tk = env[[sym[1]]]
  if (!(toupper(col.name) %in% colnames(tk)) ) {stop ('cols are not in env') }
  ind = index(tk)
  for (i in 1: len(sym)) {tmp[[i]] = env[[sym[i]]][ , toupper(col.name) ] }
  tmp = abind(tmp, along = 2)
  colnames(tmp) = sym
  tmp = make.xts(tmp, order.by = ind)
  tmp
}


bt.chart.prep <- function (x, benchmark, time.windows) {
  tmp.bind = NULL
  for (i in ls(x)) {tmp.bind = merge.xts(tmp.bind, x[[i]]$ret)}
  tmp.bind = na.omit(merge.xts(tmp.bind, benchmark)[time.windows])
  tmp.bind
}


bt.env.tranc <- function (env, frm.date) {
  bt.env.clean(env)
  tks = ls(env)
  str = paste0(frm.date, '::')
  for (i in tks) {env[[i]] = env[[i]][str]}
}

bt.env.clean <- function (env) {
  rm(dates , envir = env)
  rm(execution.price , envir = env)
  rm(prices , envir = env)
  rm(symbolnames , envir = env)
  rm(weight , envir = env)
}

bt.env.bind <- function (clean.env1, clean.env2) {
  ls1 = ls(clean.env1)
  ls2 = ls(clean.env2)
  if (!identical(ls1, ls2)) {stop ('env are not identical')}
  date.idx1 = clean.env1[[ls1[1]]] %>% index %>% as.Date
  date.idx2 = clean.env2[[ls2[1]]] %>% index %>% as.Date
  if (any(date.idx1 != date.idx2)) {stop ('date index are not identical')}
  for (i in ls1) {clean.env1[[i]] = cbind.xts(clean.env1[[i]], clean.env2[[i]])}
  return(clean.env1)
}


############### Risk tools ################
risk.ES.helper <- function (ts.data, stress.ts, sec.ts, stress.sec.ts, position, nday=c(1, 5, 20, 5), probs = c(0.05, 0.05, 0.05, 0.02), alpha.flag = F ) {
  hl = list()
  
  var.1d.95 = risk.ES(ts.data, sec.ts, position, nday[1],  probs[1], alpha.flag = alpha.flag)
  var.5d.95 = risk.ES(ts.data, sec.ts, position, nday[2],  probs[2], alpha.flag = alpha.flag)
  var.20d.95 = risk.ES(ts.data, sec.ts, position, nday[3],  probs[2], alpha.flag = alpha.flag)
  svar.5d.98 = risk.ES(stress.ts, stress.sec.ts, position, nday[4],  probs[4], alpha.flag = alpha.flag)
  
  hl$var= cbind( var.1d.95$var, var.5d.95$var, var.20d.95$var, svar.5d.98$var)
  hl$tail.loss.strip = cbind( var.1d.95$tail.strip, var.5d.95$tail.strip, var.20d.95$tail.strip, svar.5d.98$tail.strip)
  risk.ES.helper = hl
  hl
}

risk.ES <- function (ts.data, sec.ts, position, nday, probs = 0.05, alpha.flag = F) {
  d = list()
  ts.data = ts.data[, match(position$Ticker, colnames(ts.data))]
  data.na.count = transform(apply(ts.data, MARGIN = 2, function(x) sum(is.na(x))))
  na.idx = which(data.na.count > 20)
  na.sec.name = position[na.idx, ]$Sector.Code
  
  if  (len(na.idx) > 0) {
    na.und = ts.data[ , na.idx]
    na.sec =  make.xts(subset(sec.ts, select = na.sec.name),  order.by = sec.ts[ , 1])
    ts.data.tmp = ts.extend(na.und, na.sec)
    ts.data[, na.idx] = ts.data.tmp 
  }
  
  delta = position$Quantity * position$Close
  sum.delta = sum(delta)
  
  if (alpha.flag == T) {
    sum.delta = sum(delta[delta >0]);
    if (count(delta[delta >0]) == 0) {sum.delta = sum(delta[delta <0]) }
  }
  
    rt.ts = ts.data/mlag(ts.data, nday) -1
  
  if (dim(rt.ts)[2] < 2) {
    pl.strp = transform(apply(rt.ts, MARGIN = 1, function(x) x *delta))
  } else {pl.strp = transform(colSums(apply(rt.ts, MARGIN = 1, function(x) x *delta), na.rm = T))}
  
  pl.strp = make.xts(pl.strp, order.by = index(rt.ts ))
  tail = pl.strp [ pl.strp  < quantile(pl.strp , na.rm = T, probs) ]
  tail = tail[order(tail)]
  wgts = expo.w (len(tail), order = T)
  col.d = paste0('d', nday, '.', (1-probs)*100, 'dol')
  col.p = paste0('d', nday, '.', (1-probs)*100, 'per')
  d$var = cbind(round(sum(tail * wgts), digits = 0), round(sum(tail * wgts)/sum.delta *100, digits = 2))
  colnames(d$var) = c(col.d, col.p)
  d$tail.strip = cbind(tail, tail/sum.delta)
  colnames(d$tail.strip) = c(col.d, col.p)
  risk.ES  = d
  risk.ES
}



risk.ri.contr <- function(ri, ret.horizon = 5, vol.horizon = 5, top.thres = 7) {
  rt5 = t(tail(ri$rt, ret.horizon)) ; vol5 = t(tail(ri$vol, vol.horizon))
  rt.chg5d = rt5[ , 5] - rt5[ , 1] ;  rt.chg1d = rt5[ , 5] - rt5[ , 4] 
  vol.chg5d = vol5[ , 5]/vol5[, 1] -1
  full.sum = cbind('ret' = rt5[ ,5], 'vol'= vol5[, 5], rt.chg1d , rt.chg5d, vol.chg5d)
  rownames (full.sum) = rownames(rt5)
  top = rt5[c(which(rank(as.vector(rt5[ ,5]), ties.method = "max") <= top.thres ) , +
                which(rank(as.vector(rt5[ ,5]), ties.method = "max") > (dim(rt5) [1] - top.thres))),  ]
  return(list(top = top, full.sum = full.sum))
}


risk.sca <- function (x, qtn = 7, tan.win = 3, neg.slope.thres = -0.02, pos.slope.thres = 0.05) {
  xe = EMA(x,3)
  qt =  quantile(x, probs = (1:qtn/qtn), na.rm = T)
  r1.tan = (rollapply(na.omit(xe), tan.win, rel.tan2) + rollapply(na.omit(xe), tan.win, rel.tan.ols))/2
  r1.tan.quick = rollapply(na.omit(xe), tan.win, rel.tan2)
  r1.implus = z.scores(diff(r1.tan.quick, 2))
  r1.varin = merge.xts(xe, r1.tan, r1.tan.quick,  r1.implus)
  colnames(r1.varin ) = c('risk.app', 'risk.tan', 'risk.tan.quick', 'implus')
  
  r1.SE = which(x > 1.96 )
  nontrade = which(r1.varin$risk.app < qt[3])
  neg.tan = which(r1.varin$risk.tan < neg.slope.thres)
  very.neg = intersect(nontrade, neg.tan )
  TR = which(r1.varin$risk.tan > neg.slope.thres)
  SR = intersect(which(r1.varin$risk.tan.quick > pos.slope.thres),  TR)
  MM = intersect(which(r1.varin$risk.app > 0), which(r1.varin$risk.tan > pos.slope.thres))
  
  con.positve =  which(consecutive.changes(SMA(r1.varin$implus, 10)/2) > 2)
  imp.neg = which(r1.varin$implus < quantile(r1.varin$implus, probs = (1:qtn/qtn), na.rm = T)[2])
  imp.pos = which(r1.varin$implus > quantile(r1.varin$implus, probs = (1:qtn/qtn), na.rm = T)[qtn-1])
  LR = intersect(con.positve, imp.neg )
  BO.up = intersect(intersect(nontrade, imp.pos), TR)
  BO.dn = intersect(nontrade,  imp.neg)
    
  sca = r1.varin[ ,1] * NA
  sca = pnorm(r1.varin$risk.app)/0.55; sca[sca < 0.5] = 0.5 #pnorm(r1.implus)/0.5 - 1
  sca[r1.SE] = 1 +  sca[r1.SE] 
  sca[neg.tan] = 0.6  + sca[neg.tan]  #dnorm(r1.varin$risk.app[neg.tan])
  sca[nontrade] = 0.5 + sca[nontrade]
  sca[very.neg] = 0.2 +  sca[very.neg]
  sca[BO.dn ] = 0.3 +  sca[BO.dn ] 
  
  sca[MM] = -0.6 + sca[MM] 
  sca[LR] = 0.5 + sca[LR]
  sca[BO.up ] = - 0.6 + sca[BO.up] 
  
  sca [sca > 2.5] = 2.5
  sca [sca < 0.8] = 0.8
  sca[is.na(sca)] = 1
  colnames(sca) = 'VaR Scalar'
  sca
}


ts.extend <- function (ticker.ts, sector.ts) {
  sec.rt = sector.ts/mlag(sector.ts, -1) 
  ticker.rt = ticker.ts/mlag(ticker.ts, -1) 
  ts.extend = ticker.rt * NA
  for (i in 1: dim(ticker.ts)[2]) {
    tmp = sec.rt[which(is.na(ticker.rt[, i])), i]
    tmp.n = rbind.xts(ticker.rt[ which(!is.na(ticker.rt[, i])), i], tmp)
    tmp.n[len(tmp.n)] = 1
    ts.extend[ ,i] = tmp.n
  }
  ts.extend = rev(cumprod(rev(ts.extend)))
  ts.extend
}


###### Debugging tools ################
ri.plot.analysis <- function(ret, vol, n.row) {
  plot(t(ret[n.row, ]) ~ t(vol[n.row, ]))
  
}


###### memory manipulation ############

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(print(object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

################################ unsorted ###############################

op.matched <- function( x, to.match) {
  row.match = match(rownames(x), rownames(to.match))
  tmp.y = to.match * 0 
  tmp.y[row.match, ] = x 
  res = tmp.y + to.match
  res
} 


rev.sharp <- function (x) {
  ret = median(x, na.rm = T)
  vol = mad(x, na.rm = T)
  tmp = ret/vol
  tmp
}


ext.num <- function (xdata, xname, index, ref.col = 2, na.prev = TRUE) {
  n.tmp = toupper(xname)
  n.list = dimnames(xdata)[[3]]
  if (len(n.list[ n.list == n.tmp ]) == 0 ) {stop ('specified field is not in the array') }
  r.tmp = apply(xdata[ , , n.tmp], MARGIN = 2, as.numeric)
  tmp = r.tmp[!is.na(r.tmp[ ,ref.col]), ]
  if (na.prev) {tmp = apply(r.tmp, MARGIN = 2, ifna.prev)} 
  tmp = make.xts(tmp, order.by = index[which(!is.na(tmp[ ,ref.col]))])
  tmp
}



find.trd.date <- function (tar.date, trd.date) {
  tmp = as.character.Date(trd.date[which.min(abs(ymd(tar.date) - ymd(trd.date)))])
}

