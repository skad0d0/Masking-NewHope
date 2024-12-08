# i = 28 #387
k = 32
q = 12289  #3329 # 12289

for i in range(1,5000):
    prob = 1 - (i*(q*q)/2**k)
    if prob > 0 and prob < 0.02:
        print("rejection prob = ", prob, i)