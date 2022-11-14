# xray.jl

########################################################
function radon(A,angles)

#   angles = angles + pi/2;

    w,h = size(A);
    
    Nr_div2 = floor(sqrt(w*w+h*h)/2);

    Nr = Int64(Nr_div2*2+1);

    L = length(angles);

    sinogram = zeros(Nr,L);

    SIN = cos.(angles);
    COS = -sin.(angles);

#  SIN = sin.(angles);
#  COS = cos.(angles);

    R = range(-0.5,stop=0.5,length=Nr)*Nr

    RSIN = zeros(Nr,L);
    RCOS = zeros(Nr,L);

    w2 = w/2;
    h2 = h/2;

      for a=1:L
        for k=1:Nr
          RSIN[k,a]=R[k]*SIN[a];
          RCOS[k,a]=R[k]*COS[a];
        end
      end

    for a=1:L
        for k=1:Nr
          for l=1:Nr

              x = RCOS[k,a] - RSIN[l,a] + w2;
              y = RSIN[k,a] + RCOS[l,a] + h2;

              x1 = floor(x);
              x2 = ceil(x);
              y1 = floor(y);
              y2 = ceil(y);

              if !(x1 <= 0 || x2 > w || y1 <= 0 || y2 > h)

                # BILINEAR INTERPOLATION - STARTS
                # x1 y1 - btm left
                # x1 y2 - top left
                # x2 y1 - btm rght
                # x2 y2 - top rght
                f11 = A[Int64(x1),Int64(y1)];
                f21 = A[Int64(x1),Int64(y2)]; # that looks like, correct way
                f12 = A[Int64(x2),Int64(y1)];
                f22 = A[Int64(x2),Int64(y2)];
                  value = 0;
                if x2==x1 && y2==y1
                  value = f11;
                elseif x2==x1
                  value = 1/(y2-y1)*( f11*(y2-y) + f22*(y-y1) );
                elseif y2==y1
                  value = 1/(x2-x1)*( f11*(x2-x) + f22*(x-x1) );
                else
                  value = 1/(x2-x1)/(y2-y1)*( f11*(x2-x)*(y2-y) +
                                              f21*(x-x1)*(y2-y) +
                                              f12*(x2-x)*(y-y1) +
                                              f22*(x-x1)*(y-y1) );
                end
                # BILINEAR INTERPOLATION - ENDS

                sinogram[k,a] += value;

              end
          end
        end
    end
    return sinogram;
end

#########################################
function iradon(A,angles) # input angles [rad]

# angles = angles + pi/2;

  N, nAngles = size(A);

  I = zeros(N,N); # reconstruction

  x = range(-0.5,stop=0.5,length=N);

  filter = abs.(range(-1, 1, N));

  # FT domain filtering
  for t=1:length(angles)
      fhat = fftshift(fft(slice(A,:,t)));
      A[:,t] = real(ifft(ifftshift(fhat.*filter)));
  end

  XCOS = zeros(N,length(angles));
  XSIN = zeros(N,length(angles));
  for k=1:N
    for a=1:length(angles)
      XCOS[k,a]=x[k]*cos(angles[a]);
      XSIN[k,a]=x[k]*sin(angles[a]);
    end
  end

  for m=1:N
    for n=1:N
      for t=1:nAngles
        r = XCOS[m,t] + XSIN[n,t];
        index = Int64(round((r/sqrt(2.0) + 0.5)*N)); # sqrt(2) magnified
        # index = Int64(round((r + 0.5)*N)); # to keep pixel resolution
          if index>0 && index<=N
            I[m,n] += A[index,t];
          end
      end
    end
  end

  return I;
end