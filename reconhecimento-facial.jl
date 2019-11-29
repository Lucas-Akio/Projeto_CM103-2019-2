using Images
using LinearAlgebra
using Plots

# Este código é baseado no banco de faces da FEI disponível em: https://fei.edu.br/~cet/facedatabase.html
# Deve-se ter todas as imagens salvas em uma pasta nomeada Faces e este arquivo no mesmo diretório que a pasta
# lista é um vetor ordenado com o número dos indivíduos que serão usados para criar o banco de dados
# n é o número de faces que serão utilizadas para o treino (de 1 a 14)
# k é o número da face a ser testada (de 1 a 14)
# i_F é o número do indivíduo que terá sua face k testada no banco de dados formado pelos indivíduos em lista
# p é o número de colunas que serão usadas da matriz U da decomposição SVD (de 1 a n)

function reconhecimento(lista, n, k; i_F = 1, p = n)
    img = load("Faces/1-01.jpg")
    t = length(img)
    R, G, B = zeros(t,length(lista)*n), zeros(t,length(lista)*n), zeros(t,length(lista)*n)
    
    for i in 1:length(lista)
        # Criando matriz do indivíduo
        dir = "Faces/"
        for j in 1:n
            # Carregando imagem
            nº = string(k, pad = 2)
            img = load("Faces/$(lista[i])-$nº.jpg")

            # Criando vetor coluna
            r, g, b = Float64.(red.(img)), Float64.(green.(img)), Float64.(blue.(img))
            r, g, b = reshape(r, t, 1), reshape(g, t, 1), reshape(b, t, 1)

            # Atribuindo coluna a matriz indivíduo
            R[:,((i-1)*n)+j], G[:,((i-1)*n)+j], B[:,((i-1)*n)+j] = r, g, b
        end
    end
    
    # Criando o vetor F da face do indivíduo
    nº = string(k, pad = 2)
    img = load("Faces/$(i_F)-$nº.jpg")

    r, g, b = Float64.(red.(img)), Float64.(green.(img)), Float64.(blue.(img))
    Fr, Fg, Fb = reshape(r, t, 1), reshape(g, t, 1), reshape(b, t, 1)
    
    for i in 1:length(lista)
        Ur, Sr, Vr = svd(R[:,(i-1)* n + 1:i*n])
        Ug, Sg, Vg = svd(G[:,(i-1)* n + 1:i*n])
        Ub, Sb, Vb = svd(B[:,(i-1)* n + 1:i*n])

        U1 = Ur[:,1:p]
        U2 = Ug[:,1:p]
        U3 = Ub[:,1:p]

        # Criando a matriz X
        Xr, Xg, Xb = zeros(p,n), zeros(p,n), zeros(p,n)
        for j in 1:n
            Xr[:,j], Xg[:,j], Xb[:,j] = U1'*R[:,(i-1)*n + j], U2'*G[:,(i-1)*n + j], U3'*B[:,(i-1)*n + j]
        end

        # Criando o vetor x
        xr, xg, xb = U1'*Fr , U2'*Fg, U3'*Fb

        # Para cada face fⱼ, é calculado a distâncias das projeções.
        for j in 1:n
            dr, dg, db = sqrt((xr-Xr[:,j])'*(xr-Xr[:,j])), sqrt((xg-Xg[:,j])'*(xg-Xg[:,j])), sqrt((xb-Xb[:,j])'*(xb-Xb[:,j]))

            if i == 1 && j == 1
                global menor_dr, menor_dg, menor_db = dr, dg, db
                global i_E, jr, jg, jb = i, j, j, j
            else
                if dr[1,1] < menor_dr[1,1]
                    menor_dr = dr
                    jr = j
                    i_E = i
                end
                if dg[1,1] < menor_dg[1,1]
                    menor_dg = dg
                    jg = j
                    i_E = i
                end
                if db[1,1] < menor_db[1,1]
                    menor_db = db
                    jb = j
                    i_E = i
                end
            end
        end
    end
    
    nº = string(jr, pad = 2)
    melhor_face = load("Faces/$(lista[i_E])-$nº.jpg")
    return melhor_face
end
