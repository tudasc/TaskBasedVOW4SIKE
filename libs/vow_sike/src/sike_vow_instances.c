// Definition of instances to run attack on SIKE
#include "types/instance.h"
#include "sike_vow_instances.h"
#if defined(p_32_20)

// original instance
instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_32_20",
     .e = 16,
     .MEMORY_LOG_SIZE = 9,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 10.,
     .PRNG_SEED = 1337,
     .delta = 14,
     .jinv = {0xB12094B4902203E9, 0x0, 0xD4A5907EE6A3B76E, 0x3},
     .E = {
         {
             .a24 = {0xB3FDF7D5B1741A56, 0x9, 0x0, 0x0},
             .xp = {0x8F158588B29F46E1, 0xC, 0x0, 0x0},
             .xq = {0xFA19DF1C8EAD197C, 0x1, 0x0, 0x0},
             .xpq = {0x6DDBBCFC44E623FC, 0xA, 0x36C7F83B0ACF3F26, 0x2},
         },
         {
             .a24 = {0x47CCDC18F176DA0F, 0xB, 0x258EAF591E22A9D9, 0xC},
             .xp = {0x54B3CCFE768D942A, 0x10, 0x9D3AFB97F123CF9C, 0x3},
             .xq = {0x7DD78D0AD0110860, 0xA, 0xC66CF48514DE696, 0x9},
             .xpq = {0xCF0882EE9857DF67, 0x1, 0x3F9C1C2698187875, 0xC},
         }}}};

// my reproduced instance
// instance_t insts_stats[NUM_INSTS_STATS] = {
//     {.MODULUS = "p_32_20",
//      .e = 16,
//      .MEMORY_LOG_SIZE = 9,
//      .ALPHA = 2.25,
//      .BETA = 10.,
//      .GAMMA = 10.,
//      .PRNG_SEED = 1337,
//      .delta = 14,
//      .jinv = {0xB12094B4902203E9, 0x0, 0xD4A5907EE6A3B76E, 0x3},
//      .E = {
//          {
//              .a24 = {0xb3fdf7d5b1741a56, 0x0000000000000009, 000000000000000000, 000000000000000000},
//              .xp = {0xbc36cf9375ffa202, 0x0000000000000005, 0xa5afdbd64d5f0e90, 0x0000000000000004},
//              .xq = {0xc3c92bc9cff34fad, 000000000000000000, 0x62ab307b5e143f4c, 0x000000000000000f},
//              .xpq = {0x838ef133c2bc32a3, 0x0000000000000012, 0xeb70021d0e4058a9, 0x000000000000000b},
//          },
//          {
//              .a24 = {0x47ccdc18f176da0f, 0x000000000000000b, 0x258eaf591e22a9d9, 0x000000000000000c},
//              .xp = {0x9bb9d81c5f9c13de, 0x0000000000000007, 0x632fa1b1dde4db57, 0x000000000000000c},
//              .xq = {0x23958622d66cf60b, 000000000000000000, 0xd759e604f9e1a154, 000000000000000000},
//              .xpq = {0x16b5828a10dbf3f8, 0x0000000000000004, 0x8d1ff5ac8294559e, 0x0000000000000011},
//              //  .a24 = {0xefe701002877533c, 0x0000000000000008, 0x108b987751f900d2, 0x0000000000000002},
//              //  .xp = {0x68265909baa8d4f5, 0x0000000000000008, 0x3843e62428201110, 0x0000000000000002},
//              //  .xq = {0x742ae63f87512e8b, 0x000000000000000f, 0x2fd00caa1f1da27f, 0x000000000000000b},
//              //  .xpq = {0x2600b6d5fa7998d5, 0x0000000000000005, 0xf90bc604521b516d, 0x0000000000000001},

//          }}}};
#elif defined(p_36_22)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_36_22",
     .e = 18,
     .MEMORY_LOG_SIZE = 10,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 10.,
     .PRNG_SEED = 1337,
     .delta = 16,
     .jinv = {0x746F73A9CFAA13E5, 0xC55, 0x7A1E90E1166968FA, 0x124},
     .E = {
         {
             .a24 = {0xD1707E4C49EAFA66, 0x90C, 0x0, 0x0},
             .xp = {0x682B62853F71D736, 0x4EA, 0x0, 0x0},
             .xq = {0x9D038855DB13E7EC, 0xB15, 0x0, 0x0},
             .xpq = {0x47F5F56308F748CF, 0x48A, 0x629A10A84F171B70, 0xDD8},
         },
         {
             .a24 = {0x278AB12BB23B5554, 0x59B, 0xA2C752E877CD7B91, 0x72F},
             .xp = {0xAECDC850C7C72C1C, 0xE23, 0xFC28CEDB420B686E, 0x25},
             .xq = {0x8A2CDD104BA6C91D, 0x42A, 0xBD769AC24549DFB3, 0xC90},
             .xpq = {0xC67C1B31A4C13D83, 0xDB6, 0xF73D3573EFEBD873, 0x7BD},
         }}}};

// instance_t insts_stats[NUM_INSTS_STATS] = {
//     {.MODULUS = "p_36_22",
//      .e = 18,
//      .MEMORY_LOG_SIZE = 10,
//      .ALPHA = 2.25,
//      .BETA = 10.,
//      .GAMMA = 10.,
//      .PRNG_SEED = 1337,
//      .delta = 16,
//      .jinv = {0x746F73A9CFAA13E5, 0xC55, 0x7A1E90E1166968FA, 0x124},
//      .E = {
//          {
//              .a24 = {0xd1707e4c49eafa66, 0x000000000000090c, 000000000000000000, 000000000000000000},
//              .xp = {0x048b6847e0d7936f, 0x0000000000000d82, 0x28afcb8305782fdd, 0x0000000000000acf},
//              .xq = {0xa40dc48aa10381aa, 0x00000000000004e2, 0xc54cb4a6fc27a7bd, 0x00000000000001d4},
//              .xpq = {0x78176a2854ed4adb, 0x0000000000000597, 0xb77e9bd1a2e3d4bb, 0x0000000000000dcd},
//          },
//          {
//              .a24 = {0x5c8012840bb8e0c4, 0x00000000000006c8, 0x5ab633f0cb7b3258, 0x0000000000000b77},
//              .xp = {0x59e2c9dccaf797f2, 0x0000000000000344, 0xecd2351c7e3748aa, 0x00000000000004a4},
//              .xq = {0x7df49c2c377256f5, 0x00000000000004b0, 0xbe8748ef30c3917f, 0x0000000000000760},
//              .xpq = {0xf53bcf6a6075135e, 0x000000000000056b, 0x942f0a9495c02e83, 0x0000000000000a82},
//          }}}};
#elif defined(p_40_25)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_40_25",
     .e = 20,
     .MEMORY_LOG_SIZE = 11,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 10.,
     .PRNG_SEED = 1337,
     .delta = 18,
     .jinv = {0xDCDA79107493A60A, 0xD2BCB, 0xACF7AE3B13E34AE0, 0x21E721},
     .E = {
         {
             .a24 = {0xB632DF5BA1B1A3C2, 0xA69DD, 0x0, 0x0},
             .xp = {0x4DAF1F61392F5978, 0x111534, 0x0, 0x0},
             .xq = {0x2A0733AE45F486E6, 0xCCDD, 0x0, 0x0},
             .xpq = {0x7B685371A07BDD29, 0x16007B, 0x3394CBA4AC4D1669, 0x1639A9},
         },
         {
             .a24 = {0xDD49AF38B4A5BD2B, 0x4A582, 0x7DA18276D387564D, 0x1AF81},
             .xp = {0x2F7960D7976D9066, 0x295049, 0x82B4DC6C6F260F50, 0x121939},
             .xq = {0x69FAA7E927EC1D1C, 0x5C874, 0x2C1B68D7598EE25A, 0xEABCF},
             .xpq = {0xCB694EDA553E1DCB, 0x25BD86, 0x698D1589EBAA7361, 0xCE136},
         }}}};
#elif defined(p_44_27)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_44_27",
     .e = 22,
     .MEMORY_LOG_SIZE = 13,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 10.,
     .PRNG_SEED = 1337,
     .delta = 20,
     .jinv = {0xF200E1EE0FDD1CEE, 0x44FCD3F, 0xEE3E8A07A06102F8, 0x833572B},
     .E = {
         {
             .a24 = {0xAEF2B01FEC7C2B43, 0x67D0C62, 0x0, 0x0},
             .xp = {0xE7917ADAC796E5BE, 0x1B3835E, 0x0, 0x0},
             .xq = {0xB32F06E2F2F9B968, 0xA8B7D4E, 0x0, 0x0},
             .xpq = {0x349C821661D535DE, 0xE284648, 0x36ECA8D118DDDEEC, 0xE38204C},
         },
         {
             .a24 = {0x43E40DE3BE321C03, 0xDB1E23B, 0x8E1D161268DD23BB, 0x7301A90},
             .xp = {0x159511BA81BDA7E4, 0xDB9E694, 0xAD79E3569171FF09, 0x6EA4460},
             .xq = {0x80C9939E1B2AF108, 0x52EB16D, 0xEC27F863ED06650, 0xAEED21C},
             .xpq = {0xE91CB69C324419CE, 0x923924D, 0xCD4D1B2BE188798E, 0x5BCD1EF},
         }}}};
#elif defined(p_48_30)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_48_30",
     .e = 24,
     .MEMORY_LOG_SIZE = 13,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 20.,
     .PRNG_SEED = 1337,
     .delta = 22,
     .jinv = {0xCC48007CE3F92495, 0x9D5F3252, 0x6A9760FC62F09CC6, 0x3DF6A7564},
     .E = {
         {
             .a24 = {0x9B34000035D7CADC, 0x2D067C6D6, 0x0, 0x0},
             .xp = {0x81F3BA5B80DE5856, 0x4D017C1B3, 0x0, 0x0},
             .xq = {0xD688D1B7F720E1C2, 0x210E1D56, 0x0, 0x0},
             .xpq = {0xB85CDE96B813AF66, 0x2489388F4, 0x9753943F2327FA08, 0x8AEE50590},
         },
         {
             .a24 = {0xAC7E61694C2F8031, 0x8F2AF7734, 0x8FC8B3D625A83506, 0x907D1C380},
             .xp = {0x8F7C2049E1D66F3A, 0x2D9B74B5A, 0xF1923BD1731245DE, 0x421399336},
             .xq = {0x6359E499FFACD1B6, 0x719B99D4C, 0xC3CE850A228D5B35, 0x7D2C993DB},
             .xpq = {0x525D110F50579B22, 0x5EBAF47C6, 0x16B8EFAB2B273200, 0x2379FD276},
         }}}};
#elif defined(p_52_33)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_52_33",
     .e = 26,
     .MEMORY_LOG_SIZE = 15,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 20.,
     .PRNG_SEED = 1337,
     .delta = 24,
     .jinv = {0x2790EF26C30A65B9, 0x12AE60A14F1, 0x1F67789C2DDD404, 0x37F221A4DF},
     .E = {
         {
             .a24 = {0x14000000019ECA40, 0xDE4179381C, 0x0, 0x0},
             .xp = {0xAD5DA7F96CF82DE1, 0xD4B8333B0F, 0x0, 0x0},
             .xq = {0xBD938564CA85B17D, 0x3B81D64F01, 0x0, 0x0},
             .xpq = {0xD753B491F91EA1AA, 0xF0E0FCE1E0, 0xE4F1717E4FF80E0E, 0x1043D9F1AD0},
         },
         {
             .a24 = {0xA3CE80974BD4E39E, 0xD6A6CFE11D, 0xF20E9C9A5140016F, 0xB50967C97B},
             .xp = {0xB9AF1ABC9D77E33C, 0x7E93F6B3E9, 0x8CB58FA0F9FCFD01, 0xAAE038BD32},
             .xq = {0xA82409AF60D29ADA, 0x13F6D71E2E, 0x27A9A4BFC09DDD0F, 0xE9671DCA2D},
             .xpq = {0x31D1A3551D9903A2, 0x15E8C180AE, 0x91C0629C11FE8A8F, 0xDC33C005F7},
         }}}};
#elif defined(p_56_35)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_56_35",
     .e = 28,
     .MEMORY_LOG_SIZE = 17,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 20.,
     .PRNG_SEED = 1337,
     .delta = 26,
     .jinv = {0xF94DBFFDE730026C, 0x38822BFA9B2A, 0x8ED014E0DEE94AAA, 0x14C2256BC3ACF6},
     .E = {
         {
             .a24 = {0xB300000000000CEF, 0x21B78258721A32, 0x0, 0x0},
             .xp = {0x3521549B7C3BEE30, 0x22CC84CEF03388, 0x0, 0x0},
             .xq = {0x7F3B05BBB9D6C2C9, 0x271D257E8294AF, 0x0, 0x0},
             .xpq = {0x82A623DB544A58FC, 0x8FFB669917274, 0x9984390F96F0F23F, 0xF80256D756419},
         },
         {
             .a24 = {0x6AEF2342B54989A3, 0xD3DBE1E03FA63, 0x23DB2B18DFEE81FD, 0x265E708F91F504},
             .xp = {0x370EC20A6C2E44ED, 0x1C8F14BCDF6E5C, 0x6D909B6A745ED39A, 0x1632490F91C258},
             .xq = {0x13E120A2F17CD994, 0xA5529B41F0514, 0xEC58635622CF6BD0, 0x144BC06EF62BE8},
             .xpq = {0x37E64D961854B095, 0x25540509740A6E, 0xA858D8E1577DFE08, 0x3CEE5D1A95380},
         }}}};
#elif defined(p_60_38)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.MODULUS = "p_60_38",
     .e = 30,
     .MEMORY_LOG_SIZE = 19,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 20.,
     .PRNG_SEED = 1337,
     .delta = 28,
     .jinv = {0x7619C456D882E9C0, 0x1DE8B2BDA5B4ECD4, 0xA2A31DBCB4597C4C, 0x202238936874C1E0},
     .E = {
         {
             .a24 = {0x9000000000000007, 0x2C80370703624224, 0x0, 0x0},
             .xp = {0x54A48B09C50C9B43, 0x2B2179279193081E, 0x0, 0x0},
             .xq = {0x200E577A4F02D07E, 0x2329EB1D6DF9A3AA, 0x0, 0x0},
             .xpq = {0x8650ADC664E1509, 0x83BAAFDAD981FED, 0xB896004D7F6D6A3D, 0x17A1C2123E9072FC},
         },
         {
             .a24 = {0x300A9D0181270102, 0x3BF7234BAF20B0E4, 0x9FD726B5999DE44E, 0x41F1C4243BBF0CDB},
             .xp = {0x7555564276CDDCE0, 0x247599AB4F4649C2, 0x31D5F38A470EAA1B, 0x851311AEA1D1104},
             .xq = {0x267A9B8509F3DBD6, 0x146E850FF6D76AFA, 0x174D91453B8FFB66, 0x2B7FC7B463F5E5B9},
             .xpq = {0xFF8C679CA032AFFC, 0x12A486F4308AC58, 0xDA25897EF73A1E7F, 0x2E21106A8F4063D9},
         }}}};
#elif defined(p_216_137)

instance_t insts_stats[NUM_INSTS_STATS] = {
    {.e = 20,
     .MEMORY_LOG_SIZE = 10,
     .ALPHA = 2.25,
     .BETA = 10.,
     .GAMMA = 20.,
     .PRNG_SEED = 1337,
     .delta = 18,
     .jinv = {0x6BE8246542864A37, 0x11E362B039D68837, 0x9BC7A67DEEC54A2E, 0x4E6B64B48B589AA2, 0x4A904043E41FA41C, 0x61B809F34CAEF5E1, 0x44A755CAA245, 0xC739FAE5D5743D2C, 0x96F1175A33853137, 0x5FDED023029DBBE5, 0x14755F7F49BCAAA1, 0x4FD0FA75FD3CF48D, 0xC24C83D46D574891, 0x7F2CD7F534C3},
     .E = {
         {
             .a24 = {0xE858, 0x0, 0x0, 0x721FE809F8000000, 0xB00349F6AB3F59A9, 0xD264A8A8BEEE8219, 0x1D9DD4F7A5DB5, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
             .xp = {0x88AEA34159076AB3, 0xA5A5CB1A903DA770, 0x16B33C999300B977, 0xAA7A6F7A88CC0624, 0xAC4EF004B7D84E09, 0x79475FE0A608CB5E, 0x1E8C483F1256F, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
             .xq = {0x3BF57331DB10FDB2, 0x936F70C4C619C7B7, 0x3D4081EE4B1462B7, 0x10E9214AB8DC1FF4, 0x35B8DB8454A82F1C, 0x706C9F7736811583, 0xBF6E146508AA, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
             .xpq = {0xE94179134A7F2DE2, 0x6673CCEC7252249, 0xC79527DE95BA4360, 0xBAC76DFBDD10CE6C, 0x3CB9CA27612F8DF9, 0xBD318A5642CD294B, 0x185034C4A42D4, 0x5BD15E78CE28EC15, 0x98241F00AA5B8A13, 0x50FCBE39A343C4EF, 0x9C76C679526A6D11, 0xF7D3C1F044C57564, 0x3FEC6B54113419CD, 0x1FC45E1431738},
         },
         {
             .a24 = {0x605AB9CC9FE6D6C7, 0xA678D2205581290F, 0xDA9DB6845C7FF497, 0x50E7879BBD30A0D7, 0x6CA3F575A8D37E0, 0x1DB395D4FF2CF06C, 0xFF7843A8619E, 0x4B5DD216D07D18F5, 0x4E859F9552292A7, 0xA5286F0ABA23D855, 0xEA106E29AB43E8CC, 0xDDE55AE059796CB5, 0xF38271742C3F4D60, 0x216CD896D147E},
             .xp = {0x8ACB2783EDA6DC9D, 0xEB56FB79524430BE, 0x208FDFEB6C1FD415, 0xDCC2151B31FDE4B1, 0x486A4DA86559C4, 0xC5EFB5657A70D23E, 0xF3CB72FA9794, 0x64159D2D589F624F, 0x7A4608160FB91307, 0x4A6B884DD9FADEA7, 0xB6F8DD1D2550239B, 0x5EA96959DDF5842F, 0x9C7D7A752F62F5B3, 0x502789C15A1F},
             .xq = {0x47508122C6E97FD0, 0x19DE7CCF3CAA0F6, 0xC7582D8790F47019, 0x8985C3D827BB6082, 0x66357E0A90FFA2CF, 0x62FD9BE26798AC17, 0x10F698EF0D685, 0xD27C4F3371BFAA23, 0x7FEE588FFF7CBFC7, 0xA676C1D8A4DA3538, 0x34D1CF8255FC237C, 0xB43D2E5596CFCB8C, 0x3BB14A9F8E79A87D, 0x140C5E7827394},
             .xpq = {0xF14246746D3FA1A4, 0x405E0AD8A0BB8679, 0x114A7CADA7AB5622, 0xEBADA4555948A5AE, 0x9AAEEA436BECC7EA, 0xFDB141A820D06FD3, 0x20B51F1D11838, 0xE6401E60DE1457A8, 0xB9F23BE75F5345D4, 0x870F9FA720D2BE41, 0x2E83BCBABE5DD1F4, 0xFB51FF83BE75B656, 0x29E127B33C211135, 0x85F1D1056A67},
         }}}};
#endif
