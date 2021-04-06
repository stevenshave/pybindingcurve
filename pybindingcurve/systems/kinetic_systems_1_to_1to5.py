import numpy as np
from scipy.integrate import solve_ivp


def system01_p_l_kd__pl(p, l, kdpl, interval=(0, 100)):
    def ode(concs, t, kdpl):
        p, l, pl = concs
        r1 = -p * l + kdpl * pl
        dpdt = r1
        dldt = r1
        dpldt = -r1
        return [dpdt, dldt, dpldt]

    ode_result = solve_ivp(
        lambda t, y: ode(y, t, kdpl), interval, [p, l, 0.0], rtol=1e-12, ptol=1e-12
    ).y[:, -1]
    return {"p": ode_result[0], "l": ode_result[1], "pl": ode_result[2]}


def multisite_1_to_1(p, l, kd1, interval=(0, 100)):
    def ode_multisite_1_to_1(concs, t, kd1):
        p, l, p1_l = concs
        r1 = -p * l + kd1 * p1_l
        dpdt = 0.0 + r1
        dp1_ldt = 0.0 - r1
        dldt = r1
        return [dpdt, dldt, dp1_ldt]

    res = solve_ivp(
        lambda t, y: ode_multisite_1_to_1(y, t, kd1), interval, [p, l, 0.0]
    ).y[2:, -1]
    return (0 + 1 * res[0]) / l


def multisite_1_to_2(p, l, kd1, kd2, interval=(0, 100)):
    def ode_multisite_1_to_2(concs, t, kd1, kd2):
        p, l, p1_l, p2_l, p1_2_l = concs
        r1 = -p * l + kd1 * p1_l
        r2 = -p * l + kd2 * p2_l
        r3 = -p1_l * l + kd2 * p1_2_l
        r4 = -p2_l * l + kd1 * p1_2_l
        dpdt = 0.0 + r1 + r2
        dp1_ldt = 0.0 - r1 + r3
        dp2_ldt = 0.0 - r2 + r4
        dp1_2_ldt = 0.0 - r3 - r4
        dldt = r1 + r2 + r3 + r4
        return [dpdt, dldt, dp1_ldt, dp2_ldt, dp1_2_ldt]

    res = solve_ivp(
        lambda t, y: ode_multisite_1_to_2(y, t, kd1, kd2),
        interval,
        [p, l, 0.0, 0.0, 0.0],
    ).y[2:, -1]
    return (1 * (sum(res[0:2])) + 2 * res[2]) / l


def multisite_1_to_3(p, l, kd1, kd2, kd3, interval=(0, 100)):
    def ode_multisite_1_to_3(concs, t, kd1, kd2, kd3):
        p, l, p1_l, p2_l, p3_l, p1_2_l, p1_3_l, p2_3_l, p1_2_3_l = concs
        r1 = -p * l + kd1 * p1_l
        r2 = -p * l + kd2 * p2_l
        r3 = -p * l + kd3 * p3_l
        r4 = -p1_l * l + kd2 * p1_2_l
        r5 = -p2_l * l + kd1 * p1_2_l
        r6 = -p1_l * l + kd3 * p1_3_l
        r7 = -p3_l * l + kd1 * p1_3_l
        r8 = -p2_l * l + kd3 * p2_3_l
        r9 = -p3_l * l + kd2 * p2_3_l
        r10 = -p1_2_l * l + kd3 * p1_2_3_l
        r11 = -p1_3_l * l + kd2 * p1_2_3_l
        r12 = -p2_3_l * l + kd1 * p1_2_3_l
        dpdt = 0.0 + r1 + r2 + r3
        dp1_ldt = 0.0 - r1 + r4 + r6
        dp2_ldt = 0.0 - r2 + r5 + r8
        dp3_ldt = 0.0 - r3 + r7 + r9
        dp1_2_ldt = 0.0 - r4 - r5 + r10
        dp1_3_ldt = 0.0 - r6 - r7 + r11
        dp2_3_ldt = 0.0 - r8 - r9 + r12
        dp1_2_3_ldt = 0.0 - r10 - r11 - r12
        dldt = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11 + r12
        return [
            dpdt,
            dldt,
            dp1_ldt,
            dp2_ldt,
            dp3_ldt,
            dp1_2_ldt,
            dp1_3_ldt,
            dp2_3_ldt,
            dp1_2_3_ldt,
        ]

    res = solve_ivp(
        lambda t, y: ode_multisite_1_to_3(y, t, kd1, kd2, kd3),
        interval,
        [p, l, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ).y[2:, -1]
    return (1 * (sum(res[0:3])) + 2 * (sum(res[3:6])) + 3 * res[6]) / l


def multisite_1_to_4(p, l, kd1, kd2, kd3, kd4, interval=(0, 100)):
    def ode_multisite_1_to_4(concs, t, kd1, kd2, kd3, kd4):
        (
            p,
            l,
            p1_l,
            p2_l,
            p3_l,
            p4_l,
            p1_2_l,
            p1_3_l,
            p1_4_l,
            p2_3_l,
            p2_4_l,
            p3_4_l,
            p1_2_3_l,
            p1_2_4_l,
            p1_3_4_l,
            p2_3_4_l,
            p1_2_3_4_l,
        ) = concs
        r1 = -p * l + kd1 * p1_l
        r2 = -p * l + kd2 * p2_l
        r3 = -p * l + kd3 * p3_l
        r4 = -p * l + kd4 * p4_l
        r5 = -p1_l * l + kd2 * p1_2_l
        r6 = -p2_l * l + kd1 * p1_2_l
        r7 = -p1_l * l + kd3 * p1_3_l
        r8 = -p3_l * l + kd1 * p1_3_l
        r9 = -p1_l * l + kd4 * p1_4_l
        r10 = -p4_l * l + kd1 * p1_4_l
        r11 = -p2_l * l + kd3 * p2_3_l
        r12 = -p3_l * l + kd2 * p2_3_l
        r13 = -p2_l * l + kd4 * p2_4_l
        r14 = -p4_l * l + kd2 * p2_4_l
        r15 = -p3_l * l + kd4 * p3_4_l
        r16 = -p4_l * l + kd3 * p3_4_l
        r17 = -p1_2_l * l + kd3 * p1_2_3_l
        r18 = -p1_3_l * l + kd2 * p1_2_3_l
        r19 = -p2_3_l * l + kd1 * p1_2_3_l
        r20 = -p1_2_l * l + kd4 * p1_2_4_l
        r21 = -p1_4_l * l + kd2 * p1_2_4_l
        r22 = -p2_4_l * l + kd1 * p1_2_4_l
        r23 = -p1_3_l * l + kd4 * p1_3_4_l
        r24 = -p1_4_l * l + kd3 * p1_3_4_l
        r25 = -p3_4_l * l + kd1 * p1_3_4_l
        r26 = -p2_3_l * l + kd4 * p2_3_4_l
        r27 = -p2_4_l * l + kd3 * p2_3_4_l
        r28 = -p3_4_l * l + kd2 * p2_3_4_l
        r29 = -p1_2_3_l * l + kd4 * p1_2_3_4_l
        r30 = -p1_2_4_l * l + kd3 * p1_2_3_4_l
        r31 = -p1_3_4_l * l + kd2 * p1_2_3_4_l
        r32 = -p2_3_4_l * l + kd1 * p1_2_3_4_l
        dpdt = 0.0 + r1 + r2 + r3 + r4
        dp1_ldt = 0.0 - r1 + r5 + r7 + r9
        dp2_ldt = 0.0 - r2 + r6 + r11 + r13
        dp3_ldt = 0.0 - r3 + r8 + r12 + r15
        dp4_ldt = 0.0 - r4 + r10 + r14 + r16
        dp1_2_ldt = 0.0 - r5 - r6 + r17 + r20
        dp1_3_ldt = 0.0 - r7 - r8 + r18 + r23
        dp1_4_ldt = 0.0 - r9 - r10 + r21 + r24
        dp2_3_ldt = 0.0 - r11 - r12 + r19 + r26
        dp2_4_ldt = 0.0 - r13 - r14 + r22 + r27
        dp3_4_ldt = 0.0 - r15 - r16 + r25 + r28
        dp1_2_3_ldt = 0.0 - r17 - r18 - r19 + r29
        dp1_2_4_ldt = 0.0 - r20 - r21 - r22 + r30
        dp1_3_4_ldt = 0.0 - r23 - r24 - r25 + r31
        dp2_3_4_ldt = 0.0 - r26 - r27 - r28 + r32
        dp1_2_3_4_ldt = 0.0 - r29 - r30 - r31 - r32
        dldt = (
            r1
            + r2
            + r3
            + r4
            + r5
            + r6
            + r7
            + r8
            + r9
            + r10
            + r11
            + r12
            + r13
            + r14
            + r15
            + r16
            + r17
            + r18
            + r19
            + r20
            + r21
            + r22
            + r23
            + r24
            + r25
            + r26
            + r27
            + r28
            + r29
            + r30
            + r31
            + r32
        )
        return [
            dpdt,
            dldt,
            dp1_ldt,
            dp2_ldt,
            dp3_ldt,
            dp4_ldt,
            dp1_2_ldt,
            dp1_3_ldt,
            dp1_4_ldt,
            dp2_3_ldt,
            dp2_4_ldt,
            dp3_4_ldt,
            dp1_2_3_ldt,
            dp1_2_4_ldt,
            dp1_3_4_ldt,
            dp2_3_4_ldt,
            dp1_2_3_4_ldt,
        ]

    res = solve_ivp(
        lambda t, y: ode_multisite_1_to_4(y, t, kd1, kd2, kd3, kd4),
        interval,
        [
            p,
            l,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    ).y[2:, -1]
    return (
        1 * (sum(res[0:4])) + 2 * (sum(res[4:10])) + 3 * (sum(res[10:14])) + 4 * res[14]
    ) / l


def multisite_1_to_5(p, l, kd1, kd2, kd3, kd4, kd5, interval=(0, 100)):
    def ode_multisite_1_to_5(concs, t, kd1, kd2, kd3, kd4, kd5):
        (
            p,
            l,
            p1_l,
            p2_l,
            p3_l,
            p4_l,
            p5_l,
            p1_2_l,
            p1_3_l,
            p1_4_l,
            p1_5_l,
            p2_3_l,
            p2_4_l,
            p2_5_l,
            p3_4_l,
            p3_5_l,
            p4_5_l,
            p1_2_3_l,
            p1_2_4_l,
            p1_2_5_l,
            p1_3_4_l,
            p1_3_5_l,
            p1_4_5_l,
            p2_3_4_l,
            p2_3_5_l,
            p2_4_5_l,
            p3_4_5_l,
            p1_2_3_4_l,
            p1_2_3_5_l,
            p1_2_4_5_l,
            p1_3_4_5_l,
            p2_3_4_5_l,
            p1_2_3_4_5_l,
        ) = concs
        r1 = -p * l + kd1 * p1_l
        r2 = -p * l + kd2 * p2_l
        r3 = -p * l + kd3 * p3_l
        r4 = -p * l + kd4 * p4_l
        r5 = -p * l + kd5 * p5_l
        r6 = -p1_l * l + kd2 * p1_2_l
        r7 = -p2_l * l + kd1 * p1_2_l
        r8 = -p1_l * l + kd3 * p1_3_l
        r9 = -p3_l * l + kd1 * p1_3_l
        r10 = -p1_l * l + kd4 * p1_4_l
        r11 = -p4_l * l + kd1 * p1_4_l
        r12 = -p1_l * l + kd5 * p1_5_l
        r13 = -p5_l * l + kd1 * p1_5_l
        r14 = -p2_l * l + kd3 * p2_3_l
        r15 = -p3_l * l + kd2 * p2_3_l
        r16 = -p2_l * l + kd4 * p2_4_l
        r17 = -p4_l * l + kd2 * p2_4_l
        r18 = -p2_l * l + kd5 * p2_5_l
        r19 = -p5_l * l + kd2 * p2_5_l
        r20 = -p3_l * l + kd4 * p3_4_l
        r21 = -p4_l * l + kd3 * p3_4_l
        r22 = -p3_l * l + kd5 * p3_5_l
        r23 = -p5_l * l + kd3 * p3_5_l
        r24 = -p4_l * l + kd5 * p4_5_l
        r25 = -p5_l * l + kd4 * p4_5_l
        r26 = -p1_2_l * l + kd3 * p1_2_3_l
        r27 = -p1_3_l * l + kd2 * p1_2_3_l
        r28 = -p2_3_l * l + kd1 * p1_2_3_l
        r29 = -p1_2_l * l + kd4 * p1_2_4_l
        r30 = -p1_4_l * l + kd2 * p1_2_4_l
        r31 = -p2_4_l * l + kd1 * p1_2_4_l
        r32 = -p1_2_l * l + kd5 * p1_2_5_l
        r33 = -p1_5_l * l + kd2 * p1_2_5_l
        r34 = -p2_5_l * l + kd1 * p1_2_5_l
        r35 = -p1_3_l * l + kd4 * p1_3_4_l
        r36 = -p1_4_l * l + kd3 * p1_3_4_l
        r37 = -p3_4_l * l + kd1 * p1_3_4_l
        r38 = -p1_3_l * l + kd5 * p1_3_5_l
        r39 = -p1_5_l * l + kd3 * p1_3_5_l
        r40 = -p3_5_l * l + kd1 * p1_3_5_l
        r41 = -p1_4_l * l + kd5 * p1_4_5_l
        r42 = -p1_5_l * l + kd4 * p1_4_5_l
        r43 = -p4_5_l * l + kd1 * p1_4_5_l
        r44 = -p2_3_l * l + kd4 * p2_3_4_l
        r45 = -p2_4_l * l + kd3 * p2_3_4_l
        r46 = -p3_4_l * l + kd2 * p2_3_4_l
        r47 = -p2_3_l * l + kd5 * p2_3_5_l
        r48 = -p2_5_l * l + kd3 * p2_3_5_l
        r49 = -p3_5_l * l + kd2 * p2_3_5_l
        r50 = -p2_4_l * l + kd5 * p2_4_5_l
        r51 = -p2_5_l * l + kd4 * p2_4_5_l
        r52 = -p4_5_l * l + kd2 * p2_4_5_l
        r53 = -p3_4_l * l + kd5 * p3_4_5_l
        r54 = -p3_5_l * l + kd4 * p3_4_5_l
        r55 = -p4_5_l * l + kd3 * p3_4_5_l
        r56 = -p1_2_3_l * l + kd4 * p1_2_3_4_l
        r57 = -p1_2_4_l * l + kd3 * p1_2_3_4_l
        r58 = -p1_3_4_l * l + kd2 * p1_2_3_4_l
        r59 = -p2_3_4_l * l + kd1 * p1_2_3_4_l
        r60 = -p1_2_3_l * l + kd5 * p1_2_3_5_l
        r61 = -p1_2_5_l * l + kd3 * p1_2_3_5_l
        r62 = -p1_3_5_l * l + kd2 * p1_2_3_5_l
        r63 = -p2_3_5_l * l + kd1 * p1_2_3_5_l
        r64 = -p1_2_4_l * l + kd5 * p1_2_4_5_l
        r65 = -p1_2_5_l * l + kd4 * p1_2_4_5_l
        r66 = -p1_4_5_l * l + kd2 * p1_2_4_5_l
        r67 = -p2_4_5_l * l + kd1 * p1_2_4_5_l
        r68 = -p1_3_4_l * l + kd5 * p1_3_4_5_l
        r69 = -p1_3_5_l * l + kd4 * p1_3_4_5_l
        r70 = -p1_4_5_l * l + kd3 * p1_3_4_5_l
        r71 = -p3_4_5_l * l + kd1 * p1_3_4_5_l
        r72 = -p2_3_4_l * l + kd5 * p2_3_4_5_l
        r73 = -p2_3_5_l * l + kd4 * p2_3_4_5_l
        r74 = -p2_4_5_l * l + kd3 * p2_3_4_5_l
        r75 = -p3_4_5_l * l + kd2 * p2_3_4_5_l
        r76 = -p1_2_3_4_l * l + kd5 * p1_2_3_4_5_l
        r77 = -p1_2_3_5_l * l + kd4 * p1_2_3_4_5_l
        r78 = -p1_2_4_5_l * l + kd3 * p1_2_3_4_5_l
        r79 = -p1_3_4_5_l * l + kd2 * p1_2_3_4_5_l
        r80 = -p2_3_4_5_l * l + kd1 * p1_2_3_4_5_l
        dpdt = 0.0 + r1 + r2 + r3 + r4 + r5
        dp1_ldt = 0.0 - r1 + r6 + r8 + r10 + r12
        dp2_ldt = 0.0 - r2 + r7 + r14 + r16 + r18
        dp3_ldt = 0.0 - r3 + r9 + r15 + r20 + r22
        dp4_ldt = 0.0 - r4 + r11 + r17 + r21 + r24
        dp5_ldt = 0.0 - r5 + r13 + r19 + r23 + r25
        dp1_2_ldt = 0.0 - r6 - r7 + r26 + r29 + r32
        dp1_3_ldt = 0.0 - r8 - r9 + r27 + r35 + r38
        dp1_4_ldt = 0.0 - r10 - r11 + r30 + r36 + r41
        dp1_5_ldt = 0.0 - r12 - r13 + r33 + r39 + r42
        dp2_3_ldt = 0.0 - r14 - r15 + r28 + r44 + r47
        dp2_4_ldt = 0.0 - r16 - r17 + r31 + r45 + r50
        dp2_5_ldt = 0.0 - r18 - r19 + r34 + r48 + r51
        dp3_4_ldt = 0.0 - r20 - r21 + r37 + r46 + r53
        dp3_5_ldt = 0.0 - r22 - r23 + r40 + r49 + r54
        dp4_5_ldt = 0.0 - r24 - r25 + r43 + r52 + r55
        dp1_2_3_ldt = 0.0 - r26 - r27 - r28 + r56 + r60
        dp1_2_4_ldt = 0.0 - r29 - r30 - r31 + r57 + r64
        dp1_2_5_ldt = 0.0 - r32 - r33 - r34 + r61 + r65
        dp1_3_4_ldt = 0.0 - r35 - r36 - r37 + r58 + r68
        dp1_3_5_ldt = 0.0 - r38 - r39 - r40 + r62 + r69
        dp1_4_5_ldt = 0.0 - r41 - r42 - r43 + r66 + r70
        dp2_3_4_ldt = 0.0 - r44 - r45 - r46 + r59 + r72
        dp2_3_5_ldt = 0.0 - r47 - r48 - r49 + r63 + r73
        dp2_4_5_ldt = 0.0 - r50 - r51 - r52 + r67 + r74
        dp3_4_5_ldt = 0.0 - r53 - r54 - r55 + r71 + r75
        dp1_2_3_4_ldt = 0.0 - r56 - r57 - r58 - r59 + r76
        dp1_2_3_5_ldt = 0.0 - r60 - r61 - r62 - r63 + r77
        dp1_2_4_5_ldt = 0.0 - r64 - r65 - r66 - r67 + r78
        dp1_3_4_5_ldt = 0.0 - r68 - r69 - r70 - r71 + r79
        dp2_3_4_5_ldt = 0.0 - r72 - r73 - r74 - r75 + r80
        dp1_2_3_4_5_ldt = 0.0 - r76 - r77 - r78 - r79 - r80
        dldt = (
            r1
            + r2
            + r3
            + r4
            + r5
            + r6
            + r7
            + r8
            + r9
            + r10
            + r11
            + r12
            + r13
            + r14
            + r15
            + r16
            + r17
            + r18
            + r19
            + r20
            + r21
            + r22
            + r23
            + r24
            + r25
            + r26
            + r27
            + r28
            + r29
            + r30
            + r31
            + r32
            + r33
            + r34
            + r35
            + r36
            + r37
            + r38
            + r39
            + r40
            + r41
            + r42
            + r43
            + r44
            + r45
            + r46
            + r47
            + r48
            + r49
            + r50
            + r51
            + r52
            + r53
            + r54
            + r55
            + r56
            + r57
            + r58
            + r59
            + r60
            + r61
            + r62
            + r63
            + r64
            + r65
            + r66
            + r67
            + r68
            + r69
            + r70
            + r71
            + r72
            + r73
            + r74
            + r75
            + r76
            + r77
            + r78
            + r79
            + r80
        )
        return [
            dpdt,
            dldt,
            dp1_ldt,
            dp2_ldt,
            dp3_ldt,
            dp4_ldt,
            dp5_ldt,
            dp1_2_ldt,
            dp1_3_ldt,
            dp1_4_ldt,
            dp1_5_ldt,
            dp2_3_ldt,
            dp2_4_ldt,
            dp2_5_ldt,
            dp3_4_ldt,
            dp3_5_ldt,
            dp4_5_ldt,
            dp1_2_3_ldt,
            dp1_2_4_ldt,
            dp1_2_5_ldt,
            dp1_3_4_ldt,
            dp1_3_5_ldt,
            dp1_4_5_ldt,
            dp2_3_4_ldt,
            dp2_3_5_ldt,
            dp2_4_5_ldt,
            dp3_4_5_ldt,
            dp1_2_3_4_ldt,
            dp1_2_3_5_ldt,
            dp1_2_4_5_ldt,
            dp1_3_4_5_ldt,
            dp2_3_4_5_ldt,
            dp1_2_3_4_5_ldt,
        ]

    res = solve_ivp(
        lambda t, y: ode_multisite_1_to_5(y, t, kd1, kd2, kd3, kd4, kd5),
        interval,
        [
            p,
            l,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    ).y[2:, -1]
    return (
        1 * (sum(res[0:5]))
        + 2 * (sum(res[5:15]))
        + 3 * (sum(res[15:25]))
        + 4 * (sum(res[25:30]))
        + 5 * res[30]
    ) / l
