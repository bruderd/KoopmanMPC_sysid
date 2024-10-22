function dlift = jacobianLift_fourier_sparser3(in1)
%JACOBIANLIFT_FOURIER_SPARSER3
%    DLIFT = JACOBIANLIFT_FOURIER_SPARSER3(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    27-Jan-2019 13:34:10

ud1 = in1(5,:);
ud2 = in1(6,:);
ud3 = in1(7,:);
x1 = in1(1,:);
x2 = in1(2,:);
xd1 = in1(3,:);
xd2 = in1(4,:);
t2 = x1.*pi.*2.0;
t3 = x2.*pi.*2.0;
t4 = cos(t2);
t5 = sin(t3);
t6 = cos(t3);
t7 = sin(t2);
t8 = xd1.*pi.*2.0;
t9 = sin(t8);
t10 = xd2.*pi.*2.0;
t11 = sin(t10);
t12 = ud1.*pi.*2.0;
t13 = sin(t12);
t14 = ud2.*pi.*2.0;
t15 = sin(t14);
t16 = ud3.*pi.*2.0;
t17 = sin(t16);
t18 = x1.*pi.*4.0;
t19 = t4.*t6.*pi.*2.0;
t20 = t6.*t7.*pi.*2.0;
t21 = t4.*t5.*pi.*2.0;
t22 = x2.*pi.*4.0;
t23 = cos(t8);
t24 = cos(t10);
t25 = cos(t12);
t26 = cos(t14);
t27 = cos(t16);
t28 = cos(t18);
t29 = sin(t18);
t30 = sin(t22);
t31 = cos(t22);
t32 = xd1.*pi.*4.0;
t33 = sin(t32);
t34 = xd2.*pi.*4.0;
t35 = sin(t34);
t36 = ud1.*pi.*4.0;
t37 = sin(t36);
t38 = ud2.*pi.*4.0;
t39 = sin(t38);
t40 = ud3.*pi.*4.0;
t41 = sin(t40);
t42 = t4.^2;
t43 = t7.^2;
t44 = x1.*pi.*6.0;
t45 = t6.^2;
t46 = t5.^2;
t47 = t4.*t6.*t9.*pi.*2.0;
t48 = t4.*t6.*t11.*pi.*2.0;
t49 = t4.*t6.*t13.*pi.*2.0;
t50 = t4.*t6.*t15.*pi.*2.0;
t51 = t4.*t6.*t17.*pi.*2.0;
t52 = t6.*t7.*t9.*pi.*2.0;
t53 = t4.*t5.*t9.*pi.*2.0;
t54 = t6.*t7.*t11.*pi.*2.0;
t55 = t4.*t5.*t11.*pi.*2.0;
t56 = t6.*t7.*t13.*pi.*2.0;
t57 = t4.*t5.*t13.*pi.*2.0;
t58 = t6.*t7.*t15.*pi.*2.0;
t59 = t4.*t5.*t15.*pi.*2.0;
t60 = t6.*t7.*t17.*pi.*2.0;
t61 = t4.*t5.*t17.*pi.*2.0;
t62 = x2.*pi.*6.0;
t63 = t4.*t6.*t23.*pi.*2.0;
t64 = t6.*t7.*t23.*pi.*2.0;
t65 = t4.*t5.*t23.*pi.*2.0;
t66 = cos(t32);
t67 = t4.*t6.*t24.*pi.*2.0;
t68 = t6.*t7.*t24.*pi.*2.0;
t69 = t4.*t5.*t24.*pi.*2.0;
t70 = cos(t34);
t71 = t4.*t6.*t25.*pi.*2.0;
t72 = t6.*t7.*t25.*pi.*2.0;
t73 = t4.*t5.*t25.*pi.*2.0;
t74 = cos(t36);
t75 = t4.*t6.*t26.*pi.*2.0;
t76 = t6.*t7.*t26.*pi.*2.0;
t77 = t4.*t5.*t26.*pi.*2.0;
t78 = cos(t38);
t79 = t4.*t6.*t27.*pi.*2.0;
t80 = t6.*t7.*t27.*pi.*2.0;
t81 = t4.*t5.*t27.*pi.*2.0;
t82 = cos(t40);
dlift = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,t28.*pi.*4.0,t21,0.0,t4.*t9.*pi.*2.0,0.0,0.0,t4.*t11.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t42.*pi.*2.0-t43.*pi.*2.0,t5.*t7.*pi.*-2.0,t7.*t9.*pi.*-2.0,t7.*t11.*pi.*-2.0,t7.*t13.*pi.*-2.0,t7.*t15.*pi.*-2.0,t7.*t17.*pi.*-2.0,t29.*pi.*-4.0,t19,0.0,0.0,0.0,0.0,0.0,0.0,-t20,0.0,t4.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t23.*pi.*-2.0,0.0,0.0,t4.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t24.*pi.*-2.0,0.0,0.0,0.0,t4.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t25.*pi.*-2.0,0.0,0.0,0.0,0.0,t4.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t26.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t4.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,pi.*cos(t44).*6.0,t5.*t28.*pi.*4.0,t4.*t30.*pi.*2.0,0.0,t9.*t28.*pi.*4.0,t53,0.0,t4.*t33.*pi.*2.0,0.0,0.0,t11.*t28.*pi.*4.0,t55,0.0,t4.*t9.*t11.*pi.*2.0,0.0,0.0,t4.*t35.*pi.*2.0,0.0,0.0,0.0,t13.*t28.*pi.*4.0,t57,0.0,t4.*t9.*t13.*pi.*2.0,0.0,0.0,t4.*t11.*t13.*pi.*2.0,0.0,0.0,0.0,t4.*t37.*pi.*2.0,0.0,0.0,0.0,0.0,t15.*t28.*pi.*4.0,t59,0.0,t4.*t9.*t15.*pi.*2.0,0.0,0.0,t4.*t11.*t15.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t15.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t39.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t17.*t28.*pi.*4.0,t61,0.0,t4.*t9.*t17.*pi.*2.0,0.0,0.0,t4.*t11.*t17.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t41.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t28.*pi.*4.0-t7.*t29.*pi.*2.0,t5.*t42.*pi.*2.0-t5.*t43.*pi.*2.0,t7.*t30.*pi.*-2.0,t9.*t42.*pi.*2.0-t9.*t43.*pi.*2.0,t5.*t7.*t9.*pi.*-2.0,t7.*t33.*pi.*-2.0,t11.*t42.*pi.*2.0-t11.*t43.*pi.*2.0,t5.*t7.*t11.*pi.*-2.0,t7.*t9.*t11.*pi.*-2.0,t7.*t35.*pi.*-2.0,t13.*t42.*pi.*2.0-t13.*t43.*pi.*2.0,t5.*t7.*t13.*pi.*-2.0,t7.*t9.*t13.*pi.*-2.0,t7.*t11.*t13.*pi.*-2.0,t7.*t37.*pi.*-2.0,t15.*t42.*pi.*2.0-t15.*t43.*pi.*2.0,t5.*t7.*t15.*pi.*-2.0,t7.*t9.*t15.*pi.*-2.0,t7.*t11.*t15.*pi.*-2.0,t7.*t13.*t15.*pi.*-2.0,t7.*t39.*pi.*-2.0,t17.*t42.*pi.*2.0-t17.*t43.*pi.*2.0,t5.*t7.*t17.*pi.*-2.0,t7.*t9.*t17.*pi.*-2.0,t7.*t11.*t17.*pi.*-2.0,t7.*t13.*t17.*pi.*-2.0,t7.*t15.*t17.*pi.*-2.0,t7.*t41.*pi.*-2.0,t4.*t28.*pi.*2.0-t7.*t29.*pi.*4.0,t5.*t29.*pi.*-4.0,t9.*t29.*pi.*-4.0,t11.*t29.*pi.*-4.0,t13.*t29.*pi.*-4.0,t15.*t29.*pi.*-4.0,t17.*t29.*pi.*-4.0,pi.*sin(t44).*-6.0,t6.*t28.*pi.*4.0,t4.*t5.*t6.*pi.*2.0,0.0,t47,0.0,0.0,t48,0.0,0.0,0.0,t49,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,t51,0.0,0.0,0.0,0.0,0.0,0.0,t6.*t42.*pi.*2.0-t6.*t43.*pi.*2.0,t5.*t6.*t7.*pi.*-2.0,-t52,-t54,-t56,-t58,-t60,t6.*t29.*pi.*-4.0,t4.*t31.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t31.*pi.*-2.0,0.0,t23.*t28.*pi.*4.0,t65,0.0,t4.*t9.*t23.*pi.*2.0,0.0,0.0,t4.*t11.*t23.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t23.*t42.*pi.*2.0-t23.*t43.*pi.*2.0,t5.*t7.*t23.*pi.*-2.0,t7.*t9.*t23.*pi.*-2.0,t7.*t11.*t23.*pi.*-2.0,t7.*t13.*t23.*pi.*-2.0,t7.*t15.*t23.*pi.*-2.0,t7.*t17.*t23.*pi.*-2.0,t23.*t29.*pi.*-4.0,t63,0.0,0.0,0.0,0.0,0.0,0.0,-t64,0.0,t4.*t66.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t66.*pi.*-2.0,0.0,0.0,t24.*t28.*pi.*4.0,t69,0.0,t4.*t9.*t24.*pi.*2.0,0.0,0.0,t4.*t11.*t24.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t24.*t42.*pi.*2.0-t24.*t43.*pi.*2.0,t5.*t7.*t24.*pi.*-2.0,t7.*t9.*t24.*pi.*-2.0,t7.*t11.*t24.*pi.*-2.0,t7.*t13.*t24.*pi.*-2.0,t7.*t15.*t24.*pi.*-2.0,t7.*t17.*t24.*pi.*-2.0,t24.*t29.*pi.*-4.0,t67,0.0,0.0,0.0,0.0,0.0,0.0,-t68,0.0,t4.*t23.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t23.*t24.*pi.*-2.0,0.0,0.0,t4.*t70.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t70.*pi.*-2.0,0.0,0.0,0.0,t25.*t28.*pi.*4.0,t73,0.0,t4.*t9.*t25.*pi.*2.0,0.0,0.0,t4.*t11.*t25.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t25.*t42.*pi.*2.0-t25.*t43.*pi.*2.0,t5.*t7.*t25.*pi.*-2.0,t7.*t9.*t25.*pi.*-2.0,t7.*t11.*t25.*pi.*-2.0,t7.*t13.*t25.*pi.*-2.0,t7.*t15.*t25.*pi.*-2.0,t7.*t17.*t25.*pi.*-2.0,t25.*t29.*pi.*-4.0,t71,0.0,0.0,0.0,0.0,0.0,0.0,-t72,0.0,t4.*t23.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t23.*t25.*pi.*-2.0,0.0,0.0,t4.*t24.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t24.*t25.*pi.*-2.0,0.0,0.0,0.0,t4.*t74.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t74.*pi.*-2.0,0.0,0.0,0.0,0.0,t26.*t28.*pi.*4.0,t77,0.0,t4.*t9.*t26.*pi.*2.0,0.0,0.0,t4.*t11.*t26.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t26.*t42.*pi.*2.0-t26.*t43.*pi.*2.0,t5.*t7.*t26.*pi.*-2.0,t7.*t9.*t26.*pi.*-2.0,t7.*t11.*t26.*pi.*-2.0,t7.*t13.*t26.*pi.*-2.0,t7.*t15.*t26.*pi.*-2.0,t7.*t17.*t26.*pi.*-2.0,t26.*t29.*pi.*-4.0,t75,0.0,0.0,0.0,0.0,0.0,0.0,-t76,0.0,t4.*t23.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t23.*t26.*pi.*-2.0,0.0,0.0,t4.*t24.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t24.*t26.*pi.*-2.0,0.0,0.0,0.0,t4.*t25.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t25.*t26.*pi.*-2.0,0.0,0.0,0.0,0.0,t4.*t78.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t78.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t27.*t28.*pi.*4.0,t81,0.0,t4.*t9.*t27.*pi.*2.0,0.0,0.0,t4.*t11.*t27.*pi.*2.0,0.0,0.0,0.0,t4.*t13.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,t4.*t15.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t4.*t17.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t27.*t42.*pi.*2.0-t27.*t43.*pi.*2.0,t5.*t7.*t27.*pi.*-2.0,t7.*t9.*t27.*pi.*-2.0,t7.*t11.*t27.*pi.*-2.0,t7.*t13.*t27.*pi.*-2.0,t7.*t15.*t27.*pi.*-2.0,t7.*t17.*t27.*pi.*-2.0,t27.*t29.*pi.*-4.0,t79,0.0,0.0,0.0,0.0,0.0,0.0,-t80,0.0,t4.*t23.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t23.*t27.*pi.*-2.0,0.0,0.0,t4.*t24.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t24.*t27.*pi.*-2.0,0.0,0.0,0.0,t4.*t25.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t25.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,t4.*t26.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t26.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t4.*t82.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*t82.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,t20,t31.*pi.*4.0,0.0,t6.*t9.*pi.*2.0,0.0,0.0,t6.*t11.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t19,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*pi.*-2.0,t45.*pi.*2.0-t46.*pi.*2.0,t5.*t9.*pi.*-2.0,t5.*t11.*pi.*-2.0,t5.*t13.*pi.*-2.0,t5.*t15.*pi.*-2.0,t5.*t17.*pi.*-2.0,-t21,t30.*pi.*-4.0,0.0,t6.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t23.*pi.*-2.0,0.0,0.0,t6.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t24.*pi.*-2.0,0.0,0.0,0.0,t6.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t25.*pi.*-2.0,0.0,0.0,0.0,0.0,t6.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t26.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t6.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,0.0,t6.*t29.*pi.*2.0,t7.*t31.*pi.*4.0,pi.*cos(t62).*6.0,0.0,t52,t9.*t31.*pi.*4.0,0.0,t6.*t33.*pi.*2.0,0.0,0.0,t54,t11.*t31.*pi.*4.0,0.0,t6.*t9.*t11.*pi.*2.0,0.0,0.0,t6.*t35.*pi.*2.0,0.0,0.0,0.0,t56,t13.*t31.*pi.*4.0,0.0,t6.*t9.*t13.*pi.*2.0,0.0,0.0,t6.*t11.*t13.*pi.*2.0,0.0,0.0,0.0,t6.*t37.*pi.*2.0,0.0,0.0,0.0,0.0,t58,t15.*t31.*pi.*4.0,0.0,t6.*t9.*t15.*pi.*2.0,0.0,0.0,t6.*t11.*t15.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t15.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t39.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t60,t17.*t31.*pi.*4.0,0.0,t6.*t9.*t17.*pi.*2.0,0.0,0.0,t6.*t11.*t17.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t17.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t41.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t6.*t7.*pi.*2.0,t4.*t31.*pi.*4.0,0.0,t47,0.0,0.0,t48,0.0,0.0,0.0,t49,0.0,0.0,0.0,0.0,t50,0.0,0.0,0.0,0.0,0.0,t51,0.0,0.0,0.0,0.0,0.0,0.0,t6.*t28.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t29.*pi.*-2.0,t7.*t45.*pi.*2.0-t7.*t46.*pi.*2.0,t5.*t30.*pi.*-2.0+t6.*t31.*pi.*4.0,t5.*t7.*t9.*pi.*-2.0,t9.*t45.*pi.*2.0-t9.*t46.*pi.*2.0,t5.*t33.*pi.*-2.0,t5.*t7.*t11.*pi.*-2.0,t11.*t45.*pi.*2.0-t11.*t46.*pi.*2.0,t5.*t9.*t11.*pi.*-2.0,t5.*t35.*pi.*-2.0,t5.*t7.*t13.*pi.*-2.0,t13.*t45.*pi.*2.0-t13.*t46.*pi.*2.0,t5.*t9.*t13.*pi.*-2.0,t5.*t11.*t13.*pi.*-2.0,t5.*t37.*pi.*-2.0,t5.*t7.*t15.*pi.*-2.0,t15.*t45.*pi.*2.0-t15.*t46.*pi.*2.0,t5.*t9.*t15.*pi.*-2.0,t5.*t11.*t15.*pi.*-2.0,t5.*t13.*t15.*pi.*-2.0,t5.*t39.*pi.*-2.0,t5.*t7.*t17.*pi.*-2.0,t17.*t45.*pi.*2.0-t17.*t46.*pi.*2.0,t5.*t9.*t17.*pi.*-2.0,t5.*t11.*t17.*pi.*-2.0,t5.*t13.*t17.*pi.*-2.0,t5.*t15.*t17.*pi.*-2.0,t5.*t41.*pi.*-2.0,t4.*t5.*t7.*pi.*-2.0,t4.*t45.*pi.*2.0-t4.*t46.*pi.*2.0,-t53,-t55,-t57,-t59,-t61,t5.*t28.*pi.*-2.0,t7.*t30.*pi.*-4.0,t5.*t30.*pi.*-4.0+t6.*t31.*pi.*2.0,t9.*t30.*pi.*-4.0,t11.*t30.*pi.*-4.0,t13.*t30.*pi.*-4.0,t15.*t30.*pi.*-4.0,t17.*t30.*pi.*-4.0,t4.*t30.*pi.*-4.0,pi.*sin(t62).*-6.0,0.0,t64,t23.*t31.*pi.*4.0,0.0,t6.*t9.*t23.*pi.*2.0,0.0,0.0,t6.*t11.*t23.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*t23.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t63,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t23.*pi.*-2.0,t23.*t45.*pi.*2.0-t23.*t46.*pi.*2.0,t5.*t9.*t23.*pi.*-2.0,t5.*t11.*t23.*pi.*-2.0,t5.*t13.*t23.*pi.*-2.0,t5.*t15.*t23.*pi.*-2.0,t5.*t17.*t23.*pi.*-2.0,-t65,t23.*t30.*pi.*-4.0,0.0,t6.*t66.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t66.*pi.*-2.0,0.0,0.0,t68,t24.*t31.*pi.*4.0,0.0,t6.*t9.*t24.*pi.*2.0,0.0,0.0,t6.*t11.*t24.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t67,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t24.*pi.*-2.0,t24.*t45.*pi.*2.0-t24.*t46.*pi.*2.0,t5.*t9.*t24.*pi.*-2.0,t5.*t11.*t24.*pi.*-2.0,t5.*t13.*t24.*pi.*-2.0,t5.*t15.*t24.*pi.*-2.0,t5.*t17.*t24.*pi.*-2.0,-t69,t24.*t30.*pi.*-4.0,0.0,t6.*t23.*t24.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t23.*t24.*pi.*-2.0,0.0,0.0,t6.*t70.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t70.*pi.*-2.0,0.0,0.0,0.0,t72,t25.*t31.*pi.*4.0,0.0,t6.*t9.*t25.*pi.*2.0,0.0,0.0,t6.*t11.*t25.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t71,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t25.*pi.*-2.0,t25.*t45.*pi.*2.0-t25.*t46.*pi.*2.0,t5.*t9.*t25.*pi.*-2.0,t5.*t11.*t25.*pi.*-2.0,t5.*t13.*t25.*pi.*-2.0,t5.*t15.*t25.*pi.*-2.0,t5.*t17.*t25.*pi.*-2.0,-t73,t25.*t30.*pi.*-4.0,0.0,t6.*t23.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t23.*t25.*pi.*-2.0,0.0,0.0,t6.*t24.*t25.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t24.*t25.*pi.*-2.0,0.0,0.0,0.0,t6.*t74.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t74.*pi.*-2.0,0.0,0.0,0.0,0.0,t76,t26.*t31.*pi.*4.0,0.0,t6.*t9.*t26.*pi.*2.0,0.0,0.0,t6.*t11.*t26.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t75,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t26.*pi.*-2.0,t26.*t45.*pi.*2.0-t26.*t46.*pi.*2.0,t5.*t9.*t26.*pi.*-2.0,t5.*t11.*t26.*pi.*-2.0,t5.*t13.*t26.*pi.*-2.0,t5.*t15.*t26.*pi.*-2.0,t5.*t17.*t26.*pi.*-2.0,-t77,t26.*t30.*pi.*-4.0,0.0,t6.*t23.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t23.*t26.*pi.*-2.0,0.0,0.0,t6.*t24.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t24.*t26.*pi.*-2.0,0.0,0.0,0.0,t6.*t25.*t26.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t25.*t26.*pi.*-2.0,0.0,0.0,0.0,0.0,t6.*t78.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t78.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t80,t27.*t31.*pi.*4.0,0.0,t6.*t9.*t27.*pi.*2.0,0.0,0.0,t6.*t11.*t27.*pi.*2.0,0.0,0.0,0.0,t6.*t13.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,t6.*t15.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,t6.*t17.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t79,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t7.*t27.*pi.*-2.0,t27.*t45.*pi.*2.0-t27.*t46.*pi.*2.0,t5.*t9.*t27.*pi.*-2.0,t5.*t11.*t27.*pi.*-2.0,t5.*t13.*t27.*pi.*-2.0,t5.*t15.*t27.*pi.*-2.0,t5.*t17.*t27.*pi.*-2.0,-t81,t27.*t30.*pi.*-4.0,0.0,t6.*t23.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t23.*t27.*pi.*-2.0,0.0,0.0,t6.*t24.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t24.*t27.*pi.*-2.0,0.0,0.0,0.0,t6.*t25.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t25.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,t6.*t26.*t27.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t26.*t27.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0,t6.*t82.*pi.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t82.*pi.*-2.0,0.0,0.0,0.0,0.0,0.0],[687,2]);
