function rhs = nls_rhs(t, ut, dummy, k);
u = ifft(ut);
rhs = -(i/2) * (k.*k).*ut + i*fft(abs(u).*abs(u).*u);
end