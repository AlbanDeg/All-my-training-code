import math

##3d
def g_3d(J,L,S):
    if J!=0:
        gj=1+((J*(J+1)-L*(L+1)+S*(S+1))/(2*J*(J+1)))
    else:
        gj=1
    m_eff=gj*math.sqrt(J*(J+1))
    m_o=2**math.sqrt(S*(S+1))
    return f"Pour J={J}, L={L} et S={S} : gj={gj} , m_o={m_o}, m_eff={m_eff}"
##4f
def g_4f(J,L,S):
    if J!=0:
        gj=1+((J*(J+1)-L*(L+1)+S*(S+1))/(2*J*(J+1)))
    else:
        gj=1
    m_eff=gj*math.sqrt(J*(J+1))
    m_o=gj*J
    return f"Pour J={J}, L={L} et S={S} : gj={gj} , m_o={m_o}, m_eff={m_eff}"

