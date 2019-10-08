function [U,D] = QRsymm(T, tol)
    N = size(T,1);
    U = eye(n,n,'like', T);
    
    p = 1;
    q = N;
    while q ~= 1
        fl = 0;
        for i = q-1:-1:1
            if (abs(T(i+1,i)) <= tol*(abs(T(i,i)) + abs(T(i+1,i+1))))
                T(i,i+1) = 0;
                T(i+1,i) = 0;

                if fl == 0
                    q = i;
                else
                    break;
                end
            else
                fl = 1;
                p = i;
            end
        end
        
        if q > 1
            %[Qm,Rm] = qr_step(T(p:q,p:q));
			Tpq = T(p:q,p:q);
			n = size(Tpq,2);
			if n == 1
				Q = 1;
				return;
			end
			d = (Tpq(n-1,n-1) - Tpq(n,n))/2;
			
			mu = Tpq(n,n) - Tpq(n,n-1)*Tpq(n,n-1)'/(d + d/norm(d)*sqrt(d*d' + Tpq(n,n-1)*Tpq(n,n-1)'));
			
			x = Tpq(1,1) - mu;
			z = Tpq(2,1);
			Qm = eye(n,n,'like', Tpq);
			for k = 1:n-1
				if abs(z) == 0 
					c = 1;
					s = 0;
				else
					if(abs(z) > abs(x))
						tau = -x/z;
						s = 1/sqrt(1+tau*tau');
						c = s*tau;
					else
						tau = -z/x;
						c = 1/sqrt(1+tau*tau');
						s = c*tau';
					end
				end
				G = eye(n,n,'like', Tpq);
				G(k:k+1,k:k+1) = [c s; -s' c];
				Tpq = G'*Tpq*G;
				Qm = Qm*G;
				
				if k < n-1
					x = Tpq(k+1,k);
					z = Tpq(k+2,k);
				end
			end
            T(p:q,p:q) = Tpq;
            
            Q = eye(n, n, 'like', T);
            Q(p:q,p:q) = Qm;
            U = U*Q;
        end
    end
    D = T;
end
