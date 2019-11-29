using Images
using Plots
using LinearAlgebra
using JLD

# imagem é o diretório da imagem
# save = true irá salvar na pasta todas as imagens geradas
# verbose = true irá a cada 10 iterações imprimir quantos % foi concluído

function compressao(imagem; save = false, verbose = false)
    img = load(imagem)
    m,n = size(img)
    m = min(m,n)
    R, G, B, A = Float64.(red.(img)), Float64.(green.(img)), Float64.(blue.(img)), Float64.(alpha.(img))
    Erro_R, Erro_G, Erro_B, Erro_A, J = zeros(m), zeros(m), zeros(m), zeros(m), zeros(m)

    Ur, Sr, Vr = svd(R)
    Ug, Sg, Vg = svd(G)
    Ub, Sb, Vb = svd(B)
    Ua, Sa, Va = svd(A)

    J = collect(1:m)
    for j in J
        U1=Ur[:,1:j]; V1=Vr[:,1:j]; S1=Sr[1:j]
        U2=Ug[:,1:j]; V2=Vg[:,1:j]; S2=Sg[1:j]
        U3=Ub[:,1:j]; V3=Vb[:,1:j]; S3=Sb[1:j]
        U4=Ua[:,1:j]; V4=Va[:,1:j]; S4=Sa[1:j]
        R1, G1, B1, A1 = U1*diagm(S1)*V1', U2*diagm(S2)*V2', U3*diagm(S3)*V3', U4*diagm(S4)*V4'
        newimg = RGBA.(clamp01nan.(R1), clamp01nan.(G1), clamp01nan.(B1),clamp01nan.(A1))
        
        if save
            nº = string(j, pad = length(digits(m)))
            save("ex$nº.png", newimg)
        end
        
        Erro_R[j], Erro_G[j], Erro_B[j], Erro_A[j] = norm(R-R1), norm(G-G1), norm(B-B1), norm(A-A1)
        
        if j % 10 == 0 && verbose
          println("$((j/m)*100)% Concluído")
        end
    end
        
    save("Erro_R.jld", "data", Erro_R)
    save("Erro_G.jld", "data", Erro_G)
    save("Erro_B.jld", "data", Erro_B)
    save("Erro_A.jld", "data", Erro_A)
    save("Erro_índice.jld", "data", J)

    plot(J,  Erro_R, color = "red",    label = "Erro_R")
    plot!(J, Erro_G, color = "green",  label = "Erro_G")
    plot!(J, Erro_B, color = "blue",   label = "Erro_B")
    plot!(J, Erro_A, color = "purple", label = "Erro_A")
end
