import sympy as sp

# This lagging symbol depends on the cache working.
class TimeSeriesSymbol(sp.Symbol):
    def __new__(cls, basename, dt=0):
        inst = sp.Symbol.__new__(cls, TimeSeriesSymbol.symbol_name(basename,dt), commutative=False)
        inst._dt = dt
        inst._basename = basename
        return inst
    
    @staticmethod
    def symbol_name(basename, dt):
        if dt == 0:
            return r'{}_t'.format(basename)
        elif dt > 0:
            return r'{}_{{t+{}}}'.format(basename,dt)
        else:
            return r'{}_{{t-{}}}'.format(basename,-dt)
    
    def lag(self, t):
        return TimeSeriesSymbol(self._basename, self._dt+t)

class LagOperator(sp.Symbol):
    def __new__(cls, name='B'):
        inst = sp.Symbol.__new__(cls, name, commutative=False)
        return inst

def ApplyLags(expr):
    def isLagOpCompatible(expr):
        if type(expr) is LagOperator:
            return True
        if type(expr) is sp.Pow:
            if type(expr.base) is LagOperator and\
               expr.exp.is_integer:
                return True
        return False

    def getLags(expr):
        if type(expr) is LagOperator:
            return -1
        if type(expr) is sp.Pow:
            return -int(expr.exp)
        raise RuntimeError("Shouldn't attempt to get lags from non-compatible expression!")

    def ApplyLagsPass(expr):
        changed = False
        if expr.func == sp.Mul:
            # Find any matching muls
            i = 0
            while i < len(expr.args)-1:
                advance = True
                if type(expr.args[i+1]) is TimeSeriesSymbol:
                    if isLagOpCompatible(expr.args[i]):
                        num_lags = getLags(expr.args[i])
                        # Lag the time series symbol and remove the lag operators.
                        
                        new_args = []
                        for j in range(0,i):
                            new_args.append(expr.args[j])
                        new_args.append(expr.args[i+1].lag(num_lags))
                        for j in range(i+2,len(expr.args)):
                            new_args.append(expr.args[j])
                        
                        expr._args = tuple(new_args)
                        advance = False
                        changed = True
                if advance:
                    i += 1
                            
        if changed:
            return changed
        for arg in expr.args:
            if ApplyLagsPass(arg):
                return True
        return False

    while ApplyLagsPass(expr):
        expr = expr.simplify()
    return expr
