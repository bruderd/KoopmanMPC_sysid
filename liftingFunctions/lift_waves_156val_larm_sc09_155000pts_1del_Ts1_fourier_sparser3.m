function fourierBasis = stateLift_fourier_sparser3(in1)
%STATELIFT_FOURIER_SPARSER3
%    FOURIERBASIS = STATELIFT_FOURIER_SPARSER3(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    14-Dec-2018 17:49:32

ud1 = in1(5,:);
ud2 = in1(6,:);
ud3 = in1(7,:);
x1 = in1(1,:);
x2 = in1(2,:);
xd1 = in1(3,:);
xd2 = in1(4,:);
t2 = x1.*pi.*2.0;
t3 = x2.*pi.*2.0;
t4 = xd1.*pi.*2.0;
t5 = xd2.*pi.*2.0;
t6 = ud1.*pi.*2.0;
t7 = ud2.*pi.*2.0;
t8 = ud3.*pi.*2.0;
t9 = sin(t2);
t10 = sin(t3);
t11 = sin(t4);
t12 = sin(t5);
t13 = sin(t6);
t14 = sin(t7);
t15 = sin(t8);
t16 = cos(t2);
t17 = x1.*pi.*4.0;
t18 = cos(t3);
t19 = x2.*pi.*4.0;
t20 = cos(t4);
t21 = xd1.*pi.*4.0;
t22 = cos(t5);
t23 = xd2.*pi.*4.0;
t24 = cos(t6);
t25 = ud1.*pi.*4.0;
t26 = cos(t7);
t27 = ud2.*pi.*4.0;
t28 = cos(t8);
t29 = ud3.*pi.*4.0;
t30 = sin(t17);
t31 = sin(t19);
t32 = sin(t21);
t33 = sin(t23);
t34 = sin(t25);
t35 = sin(t27);
t36 = sin(t29);
t37 = cos(t17);
t38 = x1.*pi.*6.0;
t39 = cos(t19);
t40 = x2.*pi.*6.0;
t41 = cos(t21);
t42 = xd1.*pi.*6.0;
t43 = cos(t23);
t44 = xd2.*pi.*6.0;
t45 = cos(t25);
t46 = ud1.*pi.*6.0;
t47 = cos(t27);
t48 = ud2.*pi.*6.0;
t49 = cos(t29);
t50 = ud3.*pi.*6.0;
fourierBasis = [x1;x2;xd1;xd2;ud1;ud2;ud3;1.0;t9;t10;t11;t12;t13;t14;t15;t16;t18;t20;t22;t24;t26;t28;t30;t9.*t10;t31;t9.*t11;t10.*t11;t32;t9.*t12;t10.*t12;t11.*t12;t33;t9.*t13;t10.*t13;t11.*t13;t12.*t13;t34;t9.*t14;t10.*t14;t11.*t14;t12.*t14;t13.*t14;t35;t9.*t15;t10.*t15;t11.*t15;t12.*t15;t13.*t15;t14.*t15;t36;t9.*t16;t10.*t16;t11.*t16;t12.*t16;t13.*t16;t14.*t16;t15.*t16;t37;t9.*t18;t10.*t18;t11.*t18;t12.*t18;t13.*t18;t14.*t18;t15.*t18;t16.*t18;t39;t9.*t20;t10.*t20;t11.*t20;t12.*t20;t13.*t20;t14.*t20;t15.*t20;t16.*t20;t18.*t20;t41;t9.*t22;t10.*t22;t11.*t22;t12.*t22;t13.*t22;t14.*t22;t15.*t22;t16.*t22;t18.*t22;t20.*t22;t43;t9.*t24;t10.*t24;t11.*t24;t12.*t24;t13.*t24;t14.*t24;t15.*t24;t16.*t24;t18.*t24;t20.*t24;t22.*t24;t45;t9.*t26;t10.*t26;t11.*t26;t12.*t26;t13.*t26;t14.*t26;t15.*t26;t16.*t26;t18.*t26;t20.*t26;t22.*t26;t24.*t26;t47;t9.*t28;t10.*t28;t11.*t28;t12.*t28;t13.*t28;t14.*t28;t15.*t28;t16.*t28;t18.*t28;t20.*t28;t22.*t28;t24.*t28;t26.*t28;t49;sin(t38);t10.*t30;t9.*t31;sin(t40);t11.*t30;t9.*t10.*t11;t11.*t31;t9.*t32;t10.*t32;sin(t42);t12.*t30;t9.*t10.*t12;t12.*t31;t9.*t11.*t12;t10.*t11.*t12;t12.*t32;t9.*t33;t10.*t33;t11.*t33;sin(t44);t13.*t30;t9.*t10.*t13;t13.*t31;t9.*t11.*t13;t10.*t11.*t13;t13.*t32;t9.*t12.*t13;t10.*t12.*t13;t11.*t12.*t13;t13.*t33;t9.*t34;t10.*t34;t11.*t34;t12.*t34;sin(t46);t14.*t30;t9.*t10.*t14;t14.*t31;t9.*t11.*t14;t10.*t11.*t14;t14.*t32;t9.*t12.*t14;t10.*t12.*t14;t11.*t12.*t14;t14.*t33;t9.*t13.*t14;t10.*t13.*t14;t11.*t13.*t14;t12.*t13.*t14;t14.*t34;t9.*t35;t10.*t35;t11.*t35;t12.*t35;t13.*t35;sin(t48);t15.*t30;t9.*t10.*t15;t15.*t31;t9.*t11.*t15;t10.*t11.*t15;t15.*t32;t9.*t12.*t15;t10.*t12.*t15;t11.*t12.*t15;t15.*t33;t9.*t13.*t15;t10.*t13.*t15;t11.*t13.*t15;t12.*t13.*t15;t15.*t34;t9.*t14.*t15;t10.*t14.*t15;t11.*t14.*t15;t12.*t14.*t15;t13.*t14.*t15;t15.*t35;t9.*t36;t10.*t36;t11.*t36;t12.*t36;t13.*t36;t14.*t36;sin(t50);t16.*t30;t9.*t10.*t16;t16.*t31;t9.*t11.*t16;t10.*t11.*t16;t16.*t32;t9.*t12.*t16;t10.*t12.*t16;t11.*t12.*t16;t16.*t33;t9.*t13.*t16;t10.*t13.*t16;t11.*t13.*t16;t12.*t13.*t16;t16.*t34;t9.*t14.*t16;t10.*t14.*t16;t11.*t14.*t16;t12.*t14.*t16;t13.*t14.*t16;t16.*t35;t9.*t15.*t16;t10.*t15.*t16;t11.*t15.*t16;t12.*t15.*t16;t13.*t15.*t16;t14.*t15.*t16;t16.*t36;t9.*t37;t10.*t37;t11.*t37;t12.*t37;t13.*t37;t14.*t37;t15.*t37;cos(t38);t18.*t30;t9.*t10.*t18;t18.*t31;t9.*t11.*t18;t10.*t11.*t18;t18.*t32;t9.*t12.*t18;t10.*t12.*t18;t11.*t12.*t18;t18.*t33;t9.*t13.*t18;t10.*t13.*t18;t11.*t13.*t18;t12.*t13.*t18;t18.*t34;t9.*t14.*t18;t10.*t14.*t18;t11.*t14.*t18;t12.*t14.*t18;t13.*t14.*t18;t18.*t35;t9.*t15.*t18;t10.*t15.*t18;t11.*t15.*t18;t12.*t15.*t18;t13.*t15.*t18;t14.*t15.*t18;t18.*t36;t9.*t16.*t18;t10.*t16.*t18;t11.*t16.*t18;t12.*t16.*t18;t13.*t16.*t18;t14.*t16.*t18;t15.*t16.*t18;t18.*t37;t9.*t39;t10.*t39;t11.*t39;t12.*t39;t13.*t39;t14.*t39;t15.*t39;t16.*t39;cos(t40);t20.*t30;t9.*t10.*t20;t20.*t31;t9.*t11.*t20;t10.*t11.*t20;t20.*t32;t9.*t12.*t20;t10.*t12.*t20;t11.*t12.*t20;t20.*t33;t9.*t13.*t20;t10.*t13.*t20;t11.*t13.*t20;t12.*t13.*t20;t20.*t34;t9.*t14.*t20;t10.*t14.*t20;t11.*t14.*t20;t12.*t14.*t20;t13.*t14.*t20;t20.*t35;t9.*t15.*t20;t10.*t15.*t20;t11.*t15.*t20;t12.*t15.*t20;t13.*t15.*t20;t14.*t15.*t20;t20.*t36;t9.*t16.*t20;t10.*t16.*t20;t11.*t16.*t20;t12.*t16.*t20;t13.*t16.*t20;t14.*t16.*t20;t15.*t16.*t20;t20.*t37;t9.*t18.*t20;t10.*t18.*t20;t11.*t18.*t20;t12.*t18.*t20;t13.*t18.*t20;t14.*t18.*t20;t15.*t18.*t20;t16.*t18.*t20;t20.*t39;t9.*t41;t10.*t41;t11.*t41;t12.*t41;t13.*t41;t14.*t41;t15.*t41;t16.*t41;t18.*t41;cos(t42);t22.*t30;t9.*t10.*t22;t22.*t31;t9.*t11.*t22;t10.*t11.*t22;t22.*t32;t9.*t12.*t22;t10.*t12.*t22;t11.*t12.*t22;t22.*t33;t9.*t13.*t22;t10.*t13.*t22;t11.*t13.*t22;t12.*t13.*t22;t22.*t34;t9.*t14.*t22;t10.*t14.*t22;t11.*t14.*t22;t12.*t14.*t22;t13.*t14.*t22;t22.*t35;t9.*t15.*t22;t10.*t15.*t22;t11.*t15.*t22;t12.*t15.*t22;t13.*t15.*t22;t14.*t15.*t22;t22.*t36;t9.*t16.*t22;t10.*t16.*t22;t11.*t16.*t22;t12.*t16.*t22;t13.*t16.*t22;t14.*t16.*t22;t15.*t16.*t22;t22.*t37;t9.*t18.*t22;t10.*t18.*t22;t11.*t18.*t22;t12.*t18.*t22;t13.*t18.*t22;t14.*t18.*t22;t15.*t18.*t22;t16.*t18.*t22;t22.*t39;t9.*t20.*t22;t10.*t20.*t22;t11.*t20.*t22;t12.*t20.*t22;t13.*t20.*t22;t14.*t20.*t22;t15.*t20.*t22;t16.*t20.*t22;t18.*t20.*t22;t22.*t41;t9.*t43;t10.*t43;t11.*t43;t12.*t43;t13.*t43;t14.*t43;t15.*t43;t16.*t43;t18.*t43;t20.*t43;cos(t44);t24.*t30;t9.*t10.*t24;t24.*t31;t9.*t11.*t24;t10.*t11.*t24;t24.*t32;t9.*t12.*t24;t10.*t12.*t24;t11.*t12.*t24;t24.*t33;t9.*t13.*t24;t10.*t13.*t24;t11.*t13.*t24;t12.*t13.*t24;t24.*t34;t9.*t14.*t24;t10.*t14.*t24;t11.*t14.*t24;t12.*t14.*t24;t13.*t14.*t24;t24.*t35;t9.*t15.*t24;t10.*t15.*t24;t11.*t15.*t24;t12.*t15.*t24;t13.*t15.*t24;t14.*t15.*t24;t24.*t36;t9.*t16.*t24;t10.*t16.*t24;t11.*t16.*t24;t12.*t16.*t24;t13.*t16.*t24;t14.*t16.*t24;t15.*t16.*t24;t24.*t37;t9.*t18.*t24;t10.*t18.*t24;t11.*t18.*t24;t12.*t18.*t24;t13.*t18.*t24;t14.*t18.*t24;t15.*t18.*t24;t16.*t18.*t24;t24.*t39;t9.*t20.*t24;t10.*t20.*t24;t11.*t20.*t24;t12.*t20.*t24;t13.*t20.*t24;t14.*t20.*t24;t15.*t20.*t24;t16.*t20.*t24;t18.*t20.*t24;t24.*t41;t9.*t22.*t24;t10.*t22.*t24;t11.*t22.*t24;t12.*t22.*t24;t13.*t22.*t24;t14.*t22.*t24;t15.*t22.*t24;t16.*t22.*t24;t18.*t22.*t24;t20.*t22.*t24;t24.*t43;t9.*t45;t10.*t45;t11.*t45;t12.*t45;t13.*t45;t14.*t45;t15.*t45;t16.*t45;t18.*t45;t20.*t45;t22.*t45;cos(t46);t26.*t30;t9.*t10.*t26;t26.*t31;t9.*t11.*t26;t10.*t11.*t26;t26.*t32;t9.*t12.*t26;t10.*t12.*t26;t11.*t12.*t26;t26.*t33;t9.*t13.*t26;t10.*t13.*t26;t11.*t13.*t26;t12.*t13.*t26;t26.*t34;t9.*t14.*t26;t10.*t14.*t26;t11.*t14.*t26;t12.*t14.*t26;t13.*t14.*t26;t26.*t35;t9.*t15.*t26;t10.*t15.*t26;t11.*t15.*t26;t12.*t15.*t26;t13.*t15.*t26;t14.*t15.*t26;t26.*t36;t9.*t16.*t26;t10.*t16.*t26;t11.*t16.*t26;t12.*t16.*t26;t13.*t16.*t26;t14.*t16.*t26;t15.*t16.*t26;t26.*t37;t9.*t18.*t26;t10.*t18.*t26;t11.*t18.*t26;t12.*t18.*t26;t13.*t18.*t26;t14.*t18.*t26;t15.*t18.*t26;t16.*t18.*t26;t26.*t39;t9.*t20.*t26;t10.*t20.*t26;t11.*t20.*t26;t12.*t20.*t26;t13.*t20.*t26;t14.*t20.*t26;t15.*t20.*t26;t16.*t20.*t26;t18.*t20.*t26;t26.*t41;t9.*t22.*t26;t10.*t22.*t26;t11.*t22.*t26;t12.*t22.*t26;t13.*t22.*t26;t14.*t22.*t26;t15.*t22.*t26;t16.*t22.*t26;t18.*t22.*t26;t20.*t22.*t26;t26.*t43;t9.*t24.*t26;t10.*t24.*t26;t11.*t24.*t26;t12.*t24.*t26;t13.*t24.*t26;t14.*t24.*t26;t15.*t24.*t26;t16.*t24.*t26;t18.*t24.*t26;t20.*t24.*t26;t22.*t24.*t26;t26.*t45;t9.*t47;t10.*t47;t11.*t47;t12.*t47;t13.*t47;t14.*t47;t15.*t47;t16.*t47;t18.*t47;t20.*t47;t22.*t47;t24.*t47;cos(t48);t28.*t30;t9.*t10.*t28;t28.*t31;t9.*t11.*t28;t10.*t11.*t28;t28.*t32;t9.*t12.*t28;t10.*t12.*t28;t11.*t12.*t28;t28.*t33;t9.*t13.*t28;t10.*t13.*t28;t11.*t13.*t28;t12.*t13.*t28;t28.*t34;t9.*t14.*t28;t10.*t14.*t28;t11.*t14.*t28;t12.*t14.*t28;t13.*t14.*t28;t28.*t35;t9.*t15.*t28;t10.*t15.*t28;t11.*t15.*t28;t12.*t15.*t28;t13.*t15.*t28;t14.*t15.*t28;t28.*t36;t9.*t16.*t28;t10.*t16.*t28;t11.*t16.*t28;t12.*t16.*t28;t13.*t16.*t28;t14.*t16.*t28;t15.*t16.*t28;t28.*t37;t9.*t18.*t28;t10.*t18.*t28;t11.*t18.*t28;t12.*t18.*t28;t13.*t18.*t28;t14.*t18.*t28;t15.*t18.*t28;t16.*t18.*t28;t28.*t39;t9.*t20.*t28;t10.*t20.*t28;t11.*t20.*t28;t12.*t20.*t28;t13.*t20.*t28;t14.*t20.*t28;t15.*t20.*t28;t16.*t20.*t28;t18.*t20.*t28;t28.*t41;t9.*t22.*t28;t10.*t22.*t28;t11.*t22.*t28;t12.*t22.*t28;t13.*t22.*t28;t14.*t22.*t28;t15.*t22.*t28;t16.*t22.*t28;t18.*t22.*t28;t20.*t22.*t28;t28.*t43;t9.*t24.*t28;t10.*t24.*t28;t11.*t24.*t28;t12.*t24.*t28;t13.*t24.*t28;t14.*t24.*t28;t15.*t24.*t28;t16.*t24.*t28;t18.*t24.*t28;t20.*t24.*t28;t22.*t24.*t28;t28.*t45;t9.*t26.*t28;t10.*t26.*t28;t11.*t26.*t28;t12.*t26.*t28;t13.*t26.*t28;t14.*t26.*t28;t15.*t26.*t28;t16.*t26.*t28;t18.*t26.*t28;t20.*t26.*t28;t22.*t26.*t28;t24.*t26.*t28;t28.*t47;t9.*t49;t10.*t49;t11.*t49;t12.*t49;t13.*t49;t14.*t49;t15.*t49;t16.*t49;t18.*t49;t20.*t49;t22.*t49;t24.*t49;t26.*t49;cos(t50)];
