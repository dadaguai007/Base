classdef Equalizer < handle
    properties
        k1;
        k2;
        k3;
        u;
        ref;
        sps;
        u1;
        u2;
        u3;
        tap;

        equalizer;

    end

    methods
        function   obj = Equalizer()
        end

        function ffe = FFE_LMS(obj, xn, dn)
            dn = dn(:);
            L1 = length(xn);
            L2 = length(dn);
            n = round(L1/obj.sps);
            debugInfo = [];
            w = zeros(obj.k1, 1);
            x = zeros(obj.k1, 1);
            yn = zeros(n,1);
            en = zeros(n,1);
            x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
            y_d = zeros(size(yn));


            for i = 1:n - 1
                x = cat(1,x(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
                yn(i) = x.'*w ;
                if i > L2
                    if yn(i) > 2
                        y_d(i) = 3;
                    elseif yn(i) > 0
                        y_d(i) = 1;
                    elseif yn(i) > -2
                        y_d(i) = -1;
                    else
                        y_d(i) = -3;
                    end
                    en(i) = y_d(i)- yn(i);
                else
                    en(i) = dn(i)- yn(i);
                end
                w = w+obj.u*en(i)*x;
            end
            en = en(:);
            yn = yn(:);

            debugInfo.y_d = dn;
            ffe = yn;
        end

        function dfe = DFE_LMS(obj, xn, dn)
            dn = dn(:);
            L1 = length(xn);
            L2 = length(dn);
            debugInfo = [];

            k_fe = obj.k1;
            k_fb = obj.k2;

            n = round(L1/obj.sps);
            yn = zeros(n, 1);
            w1 = zeros(k_fe, 1);
            w2 = zeros(k_fb, 1);
            x1 = zeros(k_fe, 1);
            x2 = zeros(k_fb, 1);
            x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));

            for i = 1:n-1
                x1 = cat(1,x1(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
                yn(i) = x1.'*w1+ x2.'*w2;
                if i+ k_fe- 1 > L2
                    if yn(i) > 2
                        y_d(i) = 3;
                    elseif yn(i) > 0
                        y_d(i) = 1;
                    elseif yn(i) > -2
                        y_d(i) = -1;
                    else
                        y_d(i) = -3;
                    end
                    x2_new = y_d(i);
                    en(i) = y_d(i)- yn(i);
                else
                    x2_new = dn(i);
                    en(i) = dn(i)- yn(i);
                end
                w1 = w1+ obj.u*en(i)*x1;
                w2 = w2+ obj.u*en(i)*x2;
                x2 = cat(1,x2(2:end),x2_new);
            end
            w = cat(1,w1,w2);
            en = en(:);
            yn = yn(:);
            debugInfo.w1 = w1;
            dfe = yn;
        end

        function vnle2 = VNLE2_LMS(obj, xn, dn)
            dn = dn(:);
            L1 = length(xn);
            L2 = length(dn);
            n = round(L1/obj.sps);
            debugInfo = [];
            xn = xn(:);
            tap_ff = obj.k1;
            tap_fb = obj.k2;
            h1 = zeros(tap_ff,1);
            x1 = zeros(tap_ff,1);
            x2 = zeros(tap_fb,1);
            x3 = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
            % x3 = xn(ref+1:end);

            m = 1;
            for k = 1:tap_fb
                for j = k:tap_fb
                    index_mat(m,:) = [k,j];
                    m = m+1;
                end
            end
            h2 = zeros(size(index_mat,1),1);

            for i = 1:n - 2
                x1 = cat(1,x1(obj.sps+1:end),x3(obj.sps*i-obj.sps+1:1:obj.sps*i));
                x2 = x1(index_mat(:,1)).*x1(index_mat(:,2));
                yn(i,:) = x1.'*h1 + x2.'*h2;
                if i + tap_ff - 1 > L2
                    [~,y_d(i)] = quantiz(yn(i),[-2, 0, 2], [-3, -1, 1, 3]);
                    en(i,:) = y_d(i)-yn(i);
                else
                    en(i,:) = dn(i)-yn(i);
                end
                h1 = h1+obj.u*en(i)*x1;
                h2 = h2+obj.u*en(i)*x2;
                w = cat(1, h1, h2);
            end
            en = en(:);
            yn = yn(:);
            debugInfo.yn = yn;
            vnle2 = yn;
        end

        function vnle3 = VNLE3_LMS(obj, xn, dn)
            dn = dn(:);
            L1 = length(xn);
            L2 = length(dn);
            debugInfo = [];

            n = round(L1/obj.sps);

            xn = xn(:);

            yn = zeros(n, 1);
            w1 = zeros(obj.k1, 1);
            x1 = zeros(obj.k1, 1);
            x2 = zeros(obj.k2, 1);
            x3 = zeros(obj.k3, 1);
            x4 = cat(1, xn(obj.ref+1:end), dn(end-obj.ref+1:end));

            m = 1;
            for k = 1:obj.k2
                for j = k:obj.k2
                    for q = j:obj.k3
                        index_mat(m,:) = [k,j,q];
                        m = m+1;
                    end
                end
            end
            w2 = zeros(size(index_mat,1),1);
            w3 = zeros(size(index_mat,1),1);

            for i = 1:n - 1
                x1 = cat(1,x1(obj.sps+1:end),x4(obj.sps*i-obj.sps+1:1:obj.sps*i));

                x2 = x1(index_mat(:,1)).*x1(index_mat(:,2));
                x3 = x1(index_mat(:,1)).*x1(index_mat(:,2)).*x1(index_mat(:,3));

                yn(i,:) = x1.'*w1 + x2.'*w2 + x3.'*w3;
                if i > L2
                    if yn(i) > 2
                        y_d(i) = 3;
                    elseif yn(i) > 0
                        y_d(i) = 1;
                    elseif yn(i) > -2
                        y_d(i) = -1;
                    else
                        y_d(i) = -3;
                    end
                    en(i,:) = y_d(i)-yn(i);
                else
                    en(i,:) = dn(i)-yn(i);
                end
                w1 = w1+obj.u1*en(i)*x1;
                w2 = w2+obj.u2*en(i)*x2;
                w3 = w3+obj.u3*en(i)*x3;
                w = cat(1, w1, w2, w3);
            end
            en = en(:);
            yn = yn(:);
            vnle3 = yn;
        end

        function vnle2_full = VNLE2_full(obj, xn, dn)

            L1 = length(xn);
            L2 = length(dn);
            L = 2* floor(L1/obj.sps/2);

            x_linear = zeros(obj.tap, 1);
            c1 = zeros(obj.tap, 1);
            c2 = zeros(obj.tap*(obj.tap+1)/2, 1);
            yn = zeros(L, 1);

            xn = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));
            m = 1;
            for i = 1:obj.tap
                for j = i:obj.tap
                    index_mat(m, :) = [i, j];
                    m = m+1;
                end
            end
            for ii = 1:L
                x_linear = cat(1,x_linear(obj.sps+1:end),xn(obj.sps*ii-obj.sps+1:1:obj.sps*ii));
                x_nonlinear =  x_linear(index_mat(:, 1)).*x_linear(index_mat(:, 2));
                yn(ii) = x_linear.'*c1 + x_nonlinear.'*c2;
                if ii+ obj.tap- 1> L2
                    [~,y_d(ii)] = quantiz(yn(ii),[-2, 0, 2], [-3, -1, 1, 3]);
                    en(ii,:) = y_d(ii) - yn(ii);
                else
                    en(ii,:) = dn(ii) - yn(ii);
                end
                c1 = c1+ 2*obj.u1* en(ii)* x_linear;
                c2 = c2+ 2*obj.u2* en(ii)* x_nonlinear*x_linear;
            end

            yn = yn(:);
            en = en(:);

            vnle2_full = yn;
        end

        function pnle2 = PNLE2_LMS(obj, xn, dn)

            L1 = length(xn);
            L2 = length(dn);
            L = 2*floor(L1/obj.sps/2);

            x_linear = zeros(obj.k1, 1);
            c1 = zeros(obj.k1, 1);
            c2 = zeros(obj.k2, 1);
            yn = zeros(L, 1);

            xn1 = cat(1, xn(obj.ref+1:end), dn(end-obj.ref+1:end));

            for i = 1:L
                x_linear = cat(1,x_linear(obj.sps+1:end),xn1(obj.sps*i-obj.sps+1:1:obj.sps*i));
                x_nonlinear = x_linear(1:obj.k2).^2;
                yn(i,:) = x_linear.'*c1 + x_nonlinear.'*c2;

                if i + obj.k1 -1 > L2
                    if yn(i) > 2
                        y_d(i) = 3;
                    elseif yn(i) > 0
                        y_d(i) = 1;
                    elseif yn(i) > -2
                        y_d(i) = -1;
                    else
                        y_d(i) = -3;
                    end
                    en(i,:) = y_d(i)-yn(i);
                else
                    en(i,:) = dn(i)-yn(i);
                end

                c1 = c1 + en(i)*obj.u1*x_linear;
                c2 = c2 + en(i)*obj.u2*x_nonlinear;
            end
            yn = yn(:);
            en = en(:);
            pnle2 = yn;
        end

        function vnle3full = VNLE3_full(obj, xn, dn)
            L1 = length(xn);
            L2 = length(dn);
            L = 2*floor(L1/obj.sps/2);

            tap_nonlinear = floor(obj.tap/2);

            x_linear = zeros(obj.tap, 1);
            c1 = zeros(obj.tap, 1);
            c2 = zeros(obj.tap*(obj.tap+1)/2, 1);
            yn = zeros(L, 1);

            xn = cat(1,xn(obj.ref+1:end),dn(end-obj.ref+1:end));

            m = 1;
            for i = 1:obj.tap
                for j = i:obj.tap
                    index_mat_nonlinear1(m, :) = [i, j];
                    m = m + 1;
                end
            end


            n = 1;
            for ii = 1:tap_nonlinear
                for jj = ii:tap_nonlinear
                    for kk = jj:tap_nonlinear
                        index_mat_nonlinear2(n, :) = [ii, jj, kk];
                        n = n+1;
                    end
                end
            end

            c3 = zeros(size(index_mat_nonlinear2, 1), 1);

            for idx = 1:L
                x_linear = cat(1,x_linear(obj.sps+1:end),xn(obj.sps*idx-obj.sps+1:1:obj.sps*idx));
                x_nonlinear1 = x_linear(index_mat_nonlinear1(:, 1)).*x_linear(index_mat_nonlinear1(:, 2));
                x_nonlinear2 = x_linear(index_mat_nonlinear2(:, 1)).*x_linear(index_mat_nonlinear2(:, 2)).*x_linear(index_mat_nonlinear2(:,3));
                yn(idx) = x_linear.'*c1 + x_nonlinear1.'*c2 + x_nonlinear2.'*c3;
                if idx+ obj.tap- 1 > L2
                    [~,y_d(idx)] = quantiz(yn(idx),[-2, 0, 2], [-3, -1, 1, 3]);
                    en(idx,:) = y_d(idx) - yn(idx);
                else
                    en(idx,:) = dn(idx) - yn(idx);
                end
                c1 = c1 + obj.u1*en(idx)*x_linear;
                c2 = c2 + obj.u2*en(idx)*x_nonlinear1;
                c3 = c3 + obj.u3*en(idx)*x_nonlinear2;
            end
            yn = yn(:);
            en = en(:);

            vnle3full = yn;
        end

        function mse = MSE(obj, equalizer)

            switch equalizer

                case 'ffe'
                    yn = obj.FFE_LMS(obj, xn, dn);

                case 'dfe'
                    yn = obj.DFE_LMS(obj, xn, dn);

                case 'vnle2_full'
                    yn = obj.VNLE2_full(obj, xn, dn);

                case 'vnle3_full'
                    yn = obj.VNLE3_full(obj, xn, dn);
            end

            for i = 1:L
                x(i) = dn(i) - yn(i);
            end
            mse = x.^2/L;
        end

    end
end