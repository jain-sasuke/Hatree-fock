% F = eSC
% F = T + V_NE + V_EE
% S = overlap matrix


% Overlap Matrix
function S = overlap(molec)

    nb = size(molec, 1);  % Get number of basis functions efficiently

    S = zeros(nb, nb);

    for i = 1:nb
        for j = 1:nb

            npi = size(molec,2);  % Get number of primitives in basis i (cell indexing)
            npj = size(molec,2);  % Get number of primitives in basis j (cell indexing)

            for k = 1:npi
                for l = 1:npj

                    % Access properties using cell indexing
                    A_i = molec(i,k).A;
                    cof_i = molec(i,k).cof;
                    alpha_i = molec(i,k).alpha;
                    cord_i = molec(i,k).cord;

                    A_j = molec(j,l).A;
                    cof_j = molec(j,l).cof;
                    alpha_j = molec(j,l).alpha;
                    cord_j = molec(j,l).cord;

                    % Calculate overlap integral components
                    p = alpha_i + alpha_j;
                    q = alpha_i * alpha_j / p;
                    r = cord_i - cord_j;
                    r2 = dot(r, r);

                    % Update overlap matrix element
                    S(i, j) = S(i, j) +( A_i * A_j * cof_i * cof_j * exp(-q * r2) * (pi / p)^(3/2));
                end
            end
        end
    end
end



% Kinetic energy
function T = KE(molec)

    nb = size(molec, 1);  % Get number of basis functions efficiently

    T = zeros(nb, nb);

    for i = 1:nb
        for j = 1:nb

            npi = size(molec,2);  % Get number of primitives
            npj = size(molec,2);  % Get number of primitives

            for k = 1:npi
                for l = 1:npj

                    % Access properties using cell indexing
                    A_i = molec(i,k).A;
                    cof_i = molec(i,k).cof;
                    alpha_i = molec(i,k).alpha;
                    cord_i = molec(i,k).cord;

                    A_j = molec(j,l).A;
                    cof_j = molec(j,l).cof;
                    alpha_j = molec(j,l).alpha;
                    cord_j = molec(j,l).cord;

                    % Calculate overlap integral components
                    N = A_i*A_j;
                    dudv = cof_i*cof_j;
                    p = alpha_i + alpha_j;
                    po = alpha_i.*cord_i + alpha_j.*cord_j;
                    pop = po/p;
                    pr = pop - cord_j;
                    rx2 = pr(1).*pr(1);
                    ry2 = pr(2).*pr(2);
                    rz2 = pr(3).*pr(3);

                    q = alpha_i * alpha_j / p;
                    r = cord_i - cord_j;
                    r2 = dot(r, r);

                    t = exp(-q*r2)*(pi/p)^(3/2)*A_i*A_j*dudv;
                    T(i,j) = T(i,j) + 3.0 * alpha_j * t;
                    T(i,j) = T(i,j) - 2.0 * (alpha_j^2) * t *(rx2 + 0.5/p);
                    T(i,j) = T(i,j) - 2.0 * (alpha_j^2) * t *(ry2 + 0.5/p);
                    T(i,j) = T(i,j) - 2.0 * (alpha_j^2) * t *(rz2 + 0.5/p);
                end
            end
        end
    end
end





% Electron Nuclear attraction

function result = boys(x, n)
    if x == 0
        result = 1.0 / (2*n + 1);
    else
        result = gammainc(n + 0.5, x) * gamma(n + 0.5) * (1.0 / (2 * x^(n + 0.5)));
    end
end


function Ven = ENA(molec, atom_cord, Z)
    natoms = size(Z,2);
    nb = size(molec, 1);  % Get number of basis functions efficiently

    Ven = zeros(nb, nb);


    for atom = 1:natoms
        for i = 1:nb
            for j = 1:nb

                npi = size(molec,2);
                npj = size(molec,2);


                for k = 1:npi
                    for l = 1:npj

                        % Access properties using cell indexing
                        A_i = molec(i,k).A;
                        cof_i = molec(i,k).cof;
                        alpha_i = molec(i,k).alpha;
                        cord_i = molec(i,k).cord;

                        A_j = molec(j,l).A;
                        cof_j = molec(j,l).cof;
                        alpha_j = molec(j,l).alpha;
                        cord_j = molec(j,l).cord;


                        N = A_i*A_j;
                        p = alpha_i + alpha_j;
                        dudv = cof_i*cof_j;
                        q = alpha_i * alpha_j / p;
                        r = cord_i - cord_j;
                        r2 = dot(r, r);


                        po = alpha_i.*cord_i + alpha_j.*cord_j;
                        pop = po/p;
                        
                        pr = pop - atom_cord(atom,:);
                        pr2 = dot(pr,pr);

                        Ven(i,j) = Ven(i,j) -  N * dudv * exp(-q*r2)*(2*pi/p) * boys(p*pr2, 0);


                    end
                end
            end
        end
    end
end




% Nuclear Nuclear Repulsion
function Vnn = NNR(atom_cord)
    nb = 2;
    Vnn = 0;
    Rx = atom_cord(1,1) - atom_cord(2,1);
    Ry = atom_cord(1,2) - atom_cord(2,2);
    Rz = atom_cord(1,3) - atom_cord(2,3);

    Rx2 = Rx*Rx;
    Ry2 = Ry*Ry;
    Rz2 = Rz*Rz;

    R = sqrt(Rx2 + Ry2 + Rz2);
    Vnn = 1.0/R;
end


% Electron Electron Repulsion

function Vee = EER(molec)
    nb = size(molec, 1);
    Vee = zeros(nb, nb, nb, nb);

    for i = 1:nb
        for j = 1:nb
            for k = 1:nb
                for l = 1:nb
                    npi = 3;

                    for ii = 1:npi
                         for jj = 1:npi
                            for kk = 1:npi
                                for ll = 1:npi

                                N = molec(i, ii).A * molec(j, jj).A * molec(k, kk).A * molec(l, ll).A;
                                didjdkdl = molec(i, ii).cof * molec(j, jj).cof * molec(k, kk).cof * molec(l, ll).cof;
                                aij = molec(i, ii).alpha + molec(j, jj).alpha;
                                akl = molec(k, kk).alpha + molec(l, ll).alpha;

                                pij = molec(i, ii).alpha.*molec(i, ii).cord + molec(j, jj).alpha.*molec(j, jj).cord;
                                pkl = molec(k, kk).alpha.*molec(k, kk).cord + molec(l, ll).alpha.*molec(l, ll).cord;
                                paij = pij/aij;
                                pakl = pkl/akl;

                                paijpakl = paij - pakl;
                                paijpakl2 = dot(paijpakl, paijpakl);
                                den = 1.0/aij + 1.0/akl;

                                qij = molec(i, ii).alpha * molec(j, jj).alpha / aij;
                                qkl = molec(k, kk).alpha * molec(l,ll).alpha / akl;

                                rij = molec(i, ii).cord - molec(j, jj).cord;
                                rkl = molec(k, kk).cord - molec(l, ll).cord;

                                rij2 = dot(rij, rij);
                                rkl2 = dot(rkl, rkl);

                            %    f1 = 2.0 * pi * pi / (aij * akl);
                            %    f2 = sqrt(pi / (aij + akl));
                            %    f3 = exp(-qij*rij2);
                            %    f4 = exp(-qkl*rkl2);

                                Vee(i,j,k,l) = Vee(i,j,k,l) + (N * didjdkdl * 2.0 * pi * pi / (aij * akl) * sqrt(pi / (aij + akl)) * exp(-qij*rij2) * exp(-qkl*rkl2) * boys(paijpakl2/den, 0));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end




% SCF Loop





function densm = com_densm(m)
    nb = size(m,1);
    densm = zeros(nb, nb);
    occ = 2.0;
    for i = 1:nb
        for j = 1:nb
            for k = 1:occ
                C = m(i, k);
                C_d = m(j,k);
                densm(i, j) = densm(i,j) + occ * C * C_d;
            end
        end
    end
end


function G = com_G(densm, eem)
    nb = 2;
    G = zeros(nb, nb);
    for i = 1:nb
        for j = 1:nb
            for k = 1:nb
                for l = 1:nb
                    dens = densm(k,l);
                    q = eem(i, j, k, l);
                    w = eem(i, l, k, j);
                    G(i,j) = G(i, j) + dens*(q - (0.5*w));
                end
            end
        end
    end
end


function e_e = com_e_e(densm, tm, enm, G)
    Ho = tm + enm;
    e_e = 0;
    nb = size(densm,1);
    for i = 1:nb
        for j = 1:nb
            e_e = e_e + (densm(i,j) * (Ho(i,j) + (0.5*G(i,j))));
        end
    end
end

function scfcc = scfc(sm, tm, enm, eem, scfp, molec)
    conv = scfp(1,1);
    step = scfp(1,2);
    e_e = 0;
    scfcc = 0;
    nb = 2;
    densm = zeros(nb, nb);
%    sm1 = sm;
%    tm1 = tm;
%    enm1 = enm;
%    eem1 = eem;

    % start scf loop
    for i = 1:step
        eeo = e_e;

        % 2ele Vee to G
        G = com_G(densm, eem);

        % Claculate F , make S unit, eigenvalue and eigen vector
        F = tm + enm + G;
        % S unit
        smi = inv(sm);
        smis = sqrtm(smi);
        % S^(-0.5)FS^(-0.5)
        Fu1 = F * smis;
        Fu = smis * Fu1;
        l = eig(Fu);
        [vec, D] = eig(Fu);
        m = smis * vec;

        %density matrix
        densm = com_densm(m);

        %electronic energy
        e_e = com_e_e(densm, tm, enm, G);

        %convergence
        if abs(e_e - eeo) < conv
            scfcc = e_e;
        end
    end
end








% Define primitive Gaussian functions


b = 60;
dist = zeros(1,b);
dist(1) = 0.5;
for i = 2:b
    dist(i) = dist(i-1) + 0.1;  % all distances in bohr unit.
end




te = zeros(1,b);
for i = 1:b

    H1pga = pg(0.3425250914E+01,  0.1543289673E+00, [0,0,0]);
    H1pgs = pg(0.6239137298E+00,  0.5353281423E+00, [0,0,0]);
    H1pgd = pg(0.1688554040E+00,  0.4446345422E+00, [0,0,0]);

    H2pga = pg(0.3425250914E+01,  0.1543289673E+00, [(i*0.1 + 0.4),0.0,0.0]);
    H2pgs = pg(0.6239137298E+00,  0.5353281423E+00, [(i*0.1 + 0.4),0,0]);
    H2pgd = pg(0.1688554040E+00,  0.4446345422E+00, [(i*0.1 + 0.4),0,0]);

    noo = 1.0;

    H1 = [H1pga, H1pgs, H1pgd];
    H2 = [H2pga, H2pgs, H2pgd];
    molec = [H1; H2];

    Z = [1.0, 1.0];
    atom_cord = [0.0,0.0,0.0; (i*0.1 + 0.4),0.0,0.0];



    sm = overlap(molec);
    tm = KE(molec);
    enmt = ENA(molec,atom_cord, Z);
    eem = EER(molec);
    vnn = NNR(atom_cord);

    scfp = [10^(-5), 20];
    e_e = scfc(sm, tm, enmt, eem, scfp, molec);

    TE = e_e + vnn;
    te(i) = TE;
    disp(eem)
    
end

plot((dist), te);
xlabel('Distance in atomic unit');
ylabel('Energy in Hatree');
title('Energy vs Dist');
grid on;
