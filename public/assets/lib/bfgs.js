// "use strict";
// (function(undefined){
// var globalScope = (typeof global !== 'undefined') ? global : this;
// var hasModule = (typeof module !== 'undefined') && module.exports;
// if(hasModule) numeric = require('./numeric');

function _bfgs(guess,obj,maxIter,eps,opt){
  var dot = numeric.dot
  var trsp = numeric.transpose;
  var mul = numeric.mul;
  var sum = numeric.sum;
  var add = numeric.add;
  var sub = numeric.sub;
  var neg = numeric.neg;
  var solve = numeric.solve;
  var addeq = numeric.addeq;
  var subeq = numeric.subeq;
  var norm2 = numeric.norm2;
  var abs = Math.abs;
  var Obj = []; 
  var W = [];
  var opt = opt || {};
  var df = obj.df;
  var f = obj.f;
  var B = numeric.identity(guess.length);
  var x = guess.slice();
  var v = {x:x,f:f(x),df:df(x)};
  var epsMin = 1e-14;
  var c1 = opt.c1||0.25;
  var c2 = opt.c2||0.75;
  function lineSearch(v,d){
    var max = 10;
    var p0 = v.f;
    var dp0 = sum(mul(v.df,d));
    if(dp0 > 0) return; 
    if(abs(dp0) < epsMin) return v; 
    var n = opt.maxTry || 25;
    function zoom(low,high){
      if(low.f === undefined) low.f = f(low.x);
      if(high.f === undefined) high.f = f(high.x);
      while(true){
        var nv = {};
        var a = nv.a = 0.5*(low.a+high.a);
        nv.s = mul(nv.a,d);
        nv.x = add(v.x,nv.s);
        var p = nv.f = f(nv.x);
        if(!n--) return lineSearchSimple(v,d);
        if(p > p0+c1*a*dp0 || p >= low.f){
          high = nv;
        }else{
          nv.df = df(nv.x);
          var dp = nv.dp = sum(mul(nv.df,d));
          if(abs(dp) <= -c2*dp0) return nv;
          if(dp*(high.a-low.a) >= 0) high = low;
          low = nv;
        }
      }
    }
    var nv = {a:1};
    var ov = {x:v.x,a:0};
    var op;
    while(true){
      nv.s = mul(nv.a,d);
      nv.x = add(v.x,nv.s);
      var a = nv.a;
      var p = nv.f = f(nv.x);
      if(p > p0+c1*a*dp0 || (op !== undefined && p > op)) return zoom(ov,nv);
      nv.df = df(nv.x);
      var dp = nv.dp = sum(mul(nv.df,d));
      if(abs(dp) <= -c2*dp0) return nv;
      if(dp >= 0) return zoom(nv,ov);
      ov.x = nv.x;
      ov.a = nv.a;
      ov.f = nv.f;
      ov.df = nv.df;
      ov.s = nv.s;
      op = p;
      nv.a = 0.5*(a+max);
      if(!n--) return lineSearchSimple(v,d);
    }
  }
  function lineSearchSimple(v,d){
    var p0 = v.f;
    var dp0 = sum(mul(v.df,d));
    if(abs(dp0) < epsMin) return v; 
    var n = opt.maxTry || 25;
    var nv = {a:dp0>0?-1:1};
    while(true){
      nv.s = mul(nv.a,d);
      nv.x = add(v.x,nv.s);
      var a = nv.a;
      var p = nv.f = f(nv.x);
      if(p < p0+c1*a*dp0){
        nv.df = df(nv.x);
        return nv;
      }
      nv.a *= 0.5;
      if(!n--) throw "too much step, during line search";
    }  
  }

  for(var i=0;i<maxIter;i++){
    var d = neg(solve(B,v.df)); // compute direction
    if(opt.maxNormDir){
      var normDir = norm2(d);
      if(normDir > opt.maxNormDir) d = mul(d,opt.maxNormDir/normDir);
    }
    
    W.push(v.x); Obj.push(v.f);
    var og = v.df;
    var nv = lineSearch(v,d); // compute step length
    if(nv === undefined){ // SemiPositiveness lost
      B = numeric.identity(guess.length);
      continue;
    }    
    if(nv == v) return [Obj, W]; // max precision
    v = nv;
    if(norm2(v.df) < eps) return [Obj, W];
    /* start: refresh Jacobian */
    var y = sub(v.df,og);
    var vs = trsp([v.s]);
    var ts = [v.s];
    var vy = trsp([y]);
    var ty = [y];
    subeq(B,mul(1/dot(ts,dot(B,vs)),dot(B,dot(vs,dot(ts,B)))));
    addeq(B,mul(1/dot(ty,vs),dot(vy,ty)));
    /* end: refresh Jacobian */
  }
  debugger;
  return [Obj, W]

}

function runBFGS2(func, w0, alpha, totalIters) {
  function f1(w){
    return func(w)[0]
  }
  function df1(w) {
    return func(w)[1]
  }
  return _bfgs(w0, {f:f1, df: df1}, totalIters, 10e-16)
}
// 
// 
// globalScope.bfgs = bfgs;
// globalScobe.runBFGS2 = runBFGS2;
// }).call(this);
